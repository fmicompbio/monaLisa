# Calculate k-mer enrichment

Given sequences, foreground/background labels and weights, calculate the
enrichment of each k-mer in foreground compared to background. This
function is called by
[`calcBinnedKmerEnr()`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedKmerEnr.md)
for each bin if `background != "model"`.

The default type of test is `"fisher"`. Alternatively, a binomial test
can be used by `test = "binomial"`. Using Fisher's exact test has the
advantage that special cases such as zero background counts are handled
without ad-hoc adjustments to the k-mer frequencies.

For `test = "fisher"`, `fisher.test` is used with
`alternative = "greater"`, making it a one-sided test for enrichment, as
is the case with the binomial test.

## Usage

``` r
.calcKmerEnrichment(k, df, test = c("fisher", "binomial"), verbose = FALSE)
```

## Arguments

- k:

  Numeric scalar giving the length of k-mers to analyze.

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

  : the weighted number of k-mer hits in foreground sequences.

- sumBackgroundWgtWithHits:

  : the weighted number of k-mer hits in background sequences.

- totalWgtForeground:

  : the total sum of weights of foreground sequences.

- totalWgtBackground:

  : the total sum of weights of background sequences.

## Details

The function works in ZOOPS mode, which means only one or zero
occurrences of a k-mer are considered per sequence. This is helpful to
reduce the impact of simple sequence repeats occurring in few sequences.
