# Bin elements of `x`.

`bin` groups elements of `x` into bins with either a constant number of
elements per bin, a constant bin width or according to user-provided bin
boundaries.

## Usage

``` r
bin(
  x,
  binmode = c("equalN", "equalWidth", "breaks"),
  nElements = round(length(x)/5),
  nBins = NULL,
  minAbsX = NULL,
  breaks = NULL,
  ...
)
```

## Arguments

- x:

  A numerical vector with the values used for binning.

- binmode:

  The algorithm to be used for binning. Possible values are: "equalN"
  (default), "equalWidth" or "breaks" (see Details).

- nElements:

  The number of elements per bin (only for `binmode="equalN"`). The
  width of bins is adjusted accordingly.

- nBins:

  The number of bins (only for `binmode="equalWidth"`). The number of
  elements per bin will be variable.

- minAbsX:

  The minimal absolute value in `x` for elements to be binned using the
  `binmode="equalN"` or `binmode="equalWidth"` (ignored for other values
  of `binmode`). Elements with `x` values in `[-minAbsX,minAbsX]` will
  be collected in a single bin.

- breaks:

  Numerical vector with bin boundaries (only for `binmode="breaks"`).
  `breaks` has to be ordered and strictly increasing, and has to be of
  length (number of bins) + 1.

- ...:

  Further arguments to be passed to
  `cut(x, breaks, include.lowest = TRUE, ...)`, such as `labels=FALSE`.

## Value

The return value from `cut(x, ...)`, typically a factor of the same
length as `x`. Binning mode, bin boundaries and the "neutral" bin are
available from `attr(..., "binmode")`, `attr(..., "breaks")` and
`attr(..., "bin0")`. For `binmode = "breaks"`, the latter will be `NA`.

## Details

Elements are binned according to the values in `x` depending on
`binmode`:

- equalN:

  Items are grouped into a variable number of bins with `nElements`
  elements each. If `minAbsX` is not `NULL`, elements with `x`-values in
  `[-minAbsX,minAbsX]` will first be collected in a single bin before
  binning the remaining elements. The boundaries of this single bin may
  be slightly adjusted in order to respect the `nElements` elements in
  the other bins.

- equalWidth:

  Items are group into `nBins` bins with a variable number of elements
  each.

- breaks:

  Items are grouped into bins using
  `cut(x, breaks, include.lowest = TRUE)`

## See also

[`cut`](https://rdrr.io/r/base/cut.html) which is used internally.

## Examples

``` r
set.seed(1)
x <- rnorm(100)
summary(bin(x, "equalN", nElements=10))
#>    [-2.21,-1.05]   (-1.05,-0.614]  (-0.614,-0.375] (-0.375,-0.0767] 
#>               10               10               10               10 
#>  (-0.0767,0.114]    (0.114,0.377]    (0.377,0.581]    (0.581,0.771] 
#>               10               10               10               10 
#>     (0.771,1.18]       (1.18,2.4] 
#>               10               10 
summary(bin(x, "equalN", nElements=10, minAbsX=0.5))
#> Warning: Zero-bin breaks (-0.621,0.576] have been adjusted by more than 20% compared to `minAbsX` to best respect `binmode`.
#> Please use `bin(..., binmode = "breaks", breaks = X)` and optionally `setZeroBin(...)` to enforce breaks as defined by `X`.
#>  [-2.21,-1.09] (-1.09,-0.621] (-0.621,0.576]  (0.576,0.769]   (0.769,1.18] 
#>             10             10             50             10             10 
#>     (1.18,2.4] 
#>             10 
summary(bin(x, "equalWidth", nBins=5))
#>  [-2.21,-1.29] (-1.29,-0.368] (-0.368,0.555]   (0.555,1.48]     (1.48,2.4] 
#>              6             24             36             28              6 
summary(bin(x, "breaks", breaks=c(-10,-1,0,1,10)))
#> [-10,-1]   (-1,0]    (0,1]   (1,10] 
#>       11       35       39       15 
```
