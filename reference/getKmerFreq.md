# Calculate observed and expected k-mer frequencies

Given a set of sequences, calculate observed and expected k-mer
frequencies. Expected frequencies are based on a Markov model of order
`MMorder`.

## Usage

``` r
getKmerFreq(
  seqs,
  kmerLen = 5,
  MMorder = 1,
  pseudocount = 1,
  zoops = TRUE,
  strata = rep(1L, length(seqs)),
  p.adjust.method = "BH",
  includeRevComp = TRUE
)
```

## Arguments

- seqs:

  Set of sequences, either a `character` vector or a
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html).

- kmerLen:

  A `numeric` scalar giving the k-mer length.

- MMorder:

  A `numeric` scalar giving the order of the Markov model used to
  calculate the expected frequencies.

- pseudocount:

  A `numeric` scalar - will be added to the observed counts for each
  k-mer to avoid zero values.

- zoops:

  A `logical` scalar. If `TRUE` (the default), only one or zero
  occurences of a k-mer are considered per sequence.

- strata:

  A `factor` or a `numeric` scalar defining the strata of sequences. A
  separate Markov model and expected k-mer frequencies are estimated for
  the set of sequences in each stratum (level in a `strata` factor). If
  `strata` is a scalar value, it will be interpreted as the number of
  strata to split the sequences into according to their CpG
  observed-over-expected counts using `kmeans(CpGoe, centers = strata)`.

- p.adjust.method:

  A character scalar selecting the p value adjustment method (used in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)).

- includeRevComp:

  A `logical` scalar. If `TRUE` (default), count k-mer occurrences in
  both `seqs` and their reverse-complement, by concatenating `seqs` and
  their reverse-complemented versions before the counting. This is
  useful if motifs can be expected to occur on any strand (e.g. DNA
  sequences of ChIP-seq peaks). If motifs are only expected on the
  forward strand (e.g. RNA sequences of CLIP-seq peaks),
  `includeRevComp = FALSE` should be used. Note that if `strata` is a
  vector of the same length as `seqs`, each reverse-complemented
  sequence will be assigned to the same stratum as the forward sequence.

## Value

A `list` with observed and expected k-mer frequencies (`freq.obs` and
`freq.exp`, respectively), and enrichment statistics for each k-mer.

## Examples

``` r
res <- getKmerFreq(seqs = c("AAAAATT", "AAATTTT"), kmerLen = 3)
names(res)
#>  [1] "freq.obs"    "freq.exp"    "log2enr"     "sqrtDelta"   "z"          
#>  [6] "p"           "padj"        "strata"      "freq.strata" "CpGoe"      
head(res$freq.obs)
#> AAA AAC AAG AAT ACA ACC 
#>   3   0   0   4   0   0 
head(res$freq.exp)
#>       AAA       AAC       AAG       AAT       ACA       ACC 
#> 1.8295455 0.1663223 0.1663223 1.4969008 0.3326446 0.3326446 
```
