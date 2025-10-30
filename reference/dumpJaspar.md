# Dump Jaspar motifs into a HOMER motif file.

Get motifs from a Jaspar database package (e.g. `JASPAR2020`) and write
them into a HOMER-compatible motif file as positional probability
matrices.

## Usage

``` r
dumpJaspar(
  filename,
  pkg = "JASPAR2020",
  opts = list(tax_group = "vertebrates"),
  pseudocount = 1,
  relScoreCutoff = 0.8,
  verbose = FALSE
)
```

## Arguments

- filename:

  Name of the output file to be created.

- pkg:

  Name of the Jaspar package to use (default: `JASPAR2020`).

- opts:

  A list with search options used in
  [`getMatrixSet`](https://rdrr.io/pkg/TFBSTools/man/getMatrixSet-methods.html).
  By default, only vertebrate motifs are included in the output using
  `opts = list(tax_group = "vertebrates")`.

- pseudocount:

  A numerical scalar with the pseudocount to be added to each element of
  the position frequency matrix extracted from Jaspar, before its
  conversion to a position probability matrix (default: 1.0).

- relScoreCutoff:

  Currently ignored. numeric(1) in \[0,1\] that sets the default motif
  log-odds score cutof to relScoreCutoff \* maximal score for each PWM
  (default: 0.8).

- verbose:

  A logical scalar. If `TRUE`, print progress messages.

## Value

`TRUE` if successful.

## See also

[`getMatrixSet`](https://rdrr.io/pkg/TFBSTools/man/getMatrixSet-methods.html)
for details on the argument `opts`.
[`homerToPFMatrixList`](https://fmicompbio.github.io/monaLisa/reference/homerToPFMatrixList.md)
to read a file with HOMER-formatted motifs into a
[`PFMatrixList`](https://rdrr.io/pkg/TFBSTools/man/XMatrixList-class.html).

## Examples

``` r
dumpJaspar(filename = tempfile(), pkg = "JASPAR2020",
           opts = list(ID = c("MA0006.1")))
#> [1] TRUE
```
