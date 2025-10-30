# Plot diagnostics of binned sequences

Plot various diagnostics of binned sequences. Three plot types are
available:

- `length`:

  plots the distribution of sequence lengths within each bin.

- `GCfrac`:

  plots the distribution of GC fractions within each bin.

- `dinucfreq`:

  plots a heatmap of the relative frequency of each dinucleotide,
  averaged across the sequences within each bin. The values are centered
  for each dinucleotide to better highlight differences between the
  bins. The average relative frequency of each dinucleotide (across the
  bins) is indicated as well.

## Usage

``` r
plotBinDiagnostics(
  seqs,
  bins,
  aspect = c("length", "GCfrac", "dinucfreq"),
  draw_quantiles = c(0.25, 0.5, 0.75),
  ...
)
```

## Arguments

- seqs:

  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object with sequences.

- bins:

  Factor of the same length and order as `seqs`, indicating the bin for
  each sequence. Typically the return value of `bin`.

- aspect:

  The diagnostic to plot. Should be one of `"length"`, `"GCfrac"` and
  `"dinucfreq"`, to plot the distribution of sequence lengths, the
  distribution of GC fractions and the average relative dinucleotide
  frequencies across the bins.

- draw_quantiles:

  For aspect=`"length"` or `"GCfrac"`, draw vertical lines at the given
  quantiles of the density estimate. If `NULL`, no quantile lines will
  be drawn.

- ...:

  Additional argument passed to `getColsByBin`.

## Value

For aspect=`"length"` or `"GCfrac"`, returns a
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object.
For aspect=`"dinucfreq"`, returns (invisibly) a
[`Heatmap-class`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap-class.html)
object.

## Examples

``` r
seqs <- Biostrings::DNAStringSet(
  vapply(1:250, function(i) paste(sample(x = c("A", "C", "G", "T"),
                                         size = round(stats::rnorm(1, 20, 5)),
                                         replace = TRUE), collapse = ""), "")
)
bins <- factor(rep(c("a", "b", "c", "d", "e"), each = 50))
plotBinDiagnostics(seqs, bins, aspect = "length")

plotBinDiagnostics(seqs, bins, aspect = "GCfrac", draw_quantiles = NULL)

plotBinDiagnostics(seqs, bins, aspect = "dinucfreq")

```
