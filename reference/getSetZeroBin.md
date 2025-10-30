# Get and set the zero bin manually

Get and set the zero bin manually

## Usage

``` r
getZeroBin(bins)

setZeroBin(bins, zeroBin)
```

## Arguments

- bins:

  Factor, typically the return value of
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md).

- zeroBin:

  Numeric or character scalar indicating the level to use as the zero
  bin, or NA.

## Value

For `getZeroBin`, the index of the level representing the zero bin. For
`setZeroBin`, a modified factor with the zero bin set to the provided
value.

## Examples

``` r
set.seed(1)
x <- rnorm(100)
bins <- bin(x, "equalN", nElements = 10, minAbsX = 0.5)
#> Warning: Zero-bin breaks (-0.621,0.576] have been adjusted by more than 20% compared to `minAbsX` to best respect `binmode`.
#> Please use `bin(..., binmode = "breaks", breaks = X)` and optionally `setZeroBin(...)` to enforce breaks as defined by `X`.
getZeroBin(bins)
#> [1] 3
bins <- setZeroBin(bins, 2)
```
