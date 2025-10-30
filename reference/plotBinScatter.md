# Scatter plot (xy-plot) of binned elements.

Plot a scatter (xy-plot) of binned elements with binning information.

## Usage

``` r
plotBinScatter(
  x,
  y,
  b,
  cols = NULL,
  xlab = deparse(substitute(x, env = as.environment(-1))),
  ylab = deparse(substitute(y, env = as.environment(-1))),
  main = "",
  legendPosition = "right",
  legend = NULL,
  legend.cex = NULL,
  ...
)
```

## Arguments

- x:

  A numerical vector with x values.

- y:

  A numerical vector with y values (the values used for binning).

- b:

  A factor that groups elements of `x,y` into bins (typically the output
  of
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md)`(y)`).

- cols:

  `NULL` or a color vector defining the colors of points. If `NULL`, the
  colors will be computed based on `b` using
  [`getColsByBin`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md)`(b)`).

- xlab:

  Label for x-axis.

- ylab:

  Label for y-axis.

- main:

  Main title.

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
  [`getColsByBin`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md)
  (only used if `cols` is `NULL`).

## Value

The generated scatter plot as a `ggplot` object.

## See also

[`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md),
[`getColsByBin`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md),
[`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html)

## Examples

``` r
set.seed(1)
x <- rnorm(100)
y <- rnorm(100)
b <- bin(y, "equalN", nElements = 10)
plotBinScatter(x, y, b)
#> Ignoring unknown labels:
#> • fill : "Bins"

```
