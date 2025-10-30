# Get colors by bin.

Get colors for elements according to their bin. Colors are assigned to
bins forming a gradient from `col1` to `col2` in the order of
`levels{b}`. `col0` is assigned to the neutral bin (attribute `""`) if
available.

## Usage

``` r
getColsByBin(
  b,
  col1 = c("#003C30", "#01665E", "#35978F", "#80CDC1", "#C7EAE5"),
  col2 = c("#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#543005"),
  col0 = "#F5F5F5"
)
```

## Arguments

- b:

  A factor that groups elements into bins (typically the output of
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md)).

- col1:

  First color.

- col2:

  Second color.

- col0:

  Neutral color.

## Value

A character vector with colors for the elements in `b`.

## See also

[`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md).

## Examples

``` r
set.seed(1)
x <- rnorm(100)
b <- bin(x, "equalN", nElements = 10)
cols <- getColsByBin(b)
```
