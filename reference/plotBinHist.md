# Histogram of binned elements.

Plot a histogram of binned elements with binning information.

## Usage

``` r
plotBinHist(
  x,
  b,
  breaks = 6 * nlevels(b),
  xlab = deparse(substitute(x, env = as.environment(-1))),
  ylab = "Frequency",
  main = "",
  legendPosition = "right",
  legend = NULL,
  legend.cex = NULL,
  ...
)
```

## Arguments

- x:

  A numerical vector with the values used for binning.

- b:

  A factor that groups elements of `x` into bins (typically the output
  of [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md)).

- breaks:

  A `numeric` scalar controlling the histogram breaks (passed to
  `geom_hist(..., bins = breaks)`).

- xlab, ylab, main:

  `character` scalars that set the x-axis label, y-axis label and the
  main title. Use `""` to suppress the label.

- legendPosition:

  A `character` scalar. If not `"none"`, draw a legend with binning
  information. The value is used to control the legend position and will
  be passed to `theme(legend.position = legendPosition)`.

- legend:

  Deprecated (ignored). Please use `legendPosition` to control the
  drawing and position of the legend.

- legend.cex:

  Deprecated (ignored). You can use
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) to set
  legend and other graphical parameters.

- ...:

  Further arguments passed to
  [`getColsByBin`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md).

## Value

The generated histogram as a `ggplot` object.

## See also

[`getColsByBin`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md),
[`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)

## Examples

``` r
set.seed(1)
x <- rnorm(100)
b <- bin(x, "equalN", nElements = 10)
plotBinHist(x, b)

```
