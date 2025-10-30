# Load output from HOMER findMotifsGenome.pl into R

Parse HOMER output files into R data structures.

## Usage

``` r
parseHomerOutput(infiles, pseudocount.log2enr = 8, p.adjust.method = "BH")
```

## Arguments

- infiles:

  HOMER output files to be parsed.

- pseudocount.log2enr:

  A numerical scalar with the pseudocount to add to foreground and
  background counts when calculating log2 motif enrichments

- p.adjust.method:

  A character scalar selecting the p value adjustment method (used in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)).

## Value

A list of nine components (`negLog10P`, `negLog10Padj`, `pearsonResid`,
`expForegroundWgtWithHits`, `log2enr`, `sumForegroundWgtWithHits` and
`sumBackgroundWgtWithHits`), seven containing each a motif (rows) by bin
(columns) matrix with raw -log10 P values, -log10 adjusted P values, the
expected number of foreground sequences with hits, the observed number
of foreground and background sequences with hits, and motif enrichments
as Pearson residuals (`pearsonResid`) and as log2 ratios (`log2enr`),
and two containing the total foreground and background weight
(`totalWgtForeground`, `totalWgtBackground`).

## Examples

``` r
outfile <- system.file("extdata", "homer_output.txt.gz",
                       package = "monaLisa")
res <- parseHomerOutput(infiles = c(bin1 = outfile))
head(res$negLog10P)
#>                            bin1
#> MA0002.2:::RUNX1     18.8006081
#> MA0003.3:::TFAP2A     5.6458283
#> MA0004.1:::Arnt       6.9530547
#> MA0006.1:::Ahr::Arnt  5.0378160
#> MA0007.3:::Ar         0.5815203
#> MA0009.2:::T          0.5172447
```
