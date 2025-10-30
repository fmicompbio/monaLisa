# Read a HOMER motif file and create a PFMatrixList

Read motifs from a file in HOMER format and create a PFMatrixList from
them.

## Usage

``` r
homerToPFMatrixList(filename, n = 100L)
```

## Arguments

- filename:

  Name of the input file with HOMER-formatted motifs.

- n:

  The number of observations (multiplied with base frequencies to create
  the number of observed bases at each position).

## Value

A
[`PFMatrixList`](https://rdrr.io/pkg/TFBSTools/man/XMatrixList-class.html)
with motifs from the file.

## See also

[`dumpJaspar`](https://fmicompbio.github.io/monaLisa/reference/dumpJaspar.md)
for writing motifs from a Jaspar database package into a file in HOMER
format.

## Examples

``` r
library(JASPAR2020)
optsL <- list(ID = c("MA0006.1"))
pfm1 <- TFBSTools::getMatrixSet(JASPAR2020, opts = optsL)
TFBSTools::Matrix(pfm1)
#> $MA0006.1
#>   [,1] [,2] [,3] [,4] [,5] [,6]
#> A    3    0    0    0    0    0
#> C    8    0   23    0    0    0
#> G    2   23    0   23    0   24
#> T   11    1    1    1   24    0
#> 

tmpfn <- tempfile()
dumpJaspar(filename = tmpfn, pkg = "JASPAR2020", opts = optsL)
#> [1] TRUE
pfm2 <- homerToPFMatrixList(tmpfn)
TFBSTools::Matrix(pfm2)
#> [[1]]
#>   [,1] [,2] [,3] [,4] [,5] [,6]
#> A   14    4    4    4    4    4
#> C   32    4   86    4    4    4
#> G   11   86    4   86    4   89
#> T   43    7    7    7   89    4
#> 

unlink(tmpfn)
```
