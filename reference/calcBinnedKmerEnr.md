# Calculate k-mer enrichment in bins of sequences.

Given a set of sequences and corresponding bins, identify enriched
k-mers (n-grams) in each bin. The sequences can be given either directly
or as genomic coordinates.

## Usage

``` r
calcBinnedKmerEnr(
  seqs,
  bins = NULL,
  kmerLen = 5,
  background = c("otherBins", "allBins", "zeroBin", "genome", "model"),
  MMorder = 1,
  test = c("fisher", "binomial"),
  includeRevComp = TRUE,
  maxFracN = 0.7,
  maxKmerSize = 3L,
  GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8),
  pseudocount.kmers = 1,
  pseudocount.log2enr = 8,
  p.adjust.method = "BH",
  genome = NULL,
  genome.regions = NULL,
  genome.oversample = 2,
  BPPARAM = SerialParam(),
  verbose = FALSE
)
```

## Arguments

- seqs:

  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object with sequences to test

- bins:

  Factor of the same length and order as `seqs`, indicating the bin for
  each sequence. Typically the return value of
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md). For
  `background = "genome"` or `background = "model"`, `bins` can be
  omitted.

- kmerLen:

  A `numeric` scalar giving the k-mer length.

- background:

  A `character` scalar specifying the background sequences to use. One
  of `"otherBins"` (default), `"allBins"`, `"zeroBin"`, `"genome"` or
  `"model"` (see "Details").

- MMorder:

  A `numeric` scalar giving the order of the Markov model used to
  calculate the expected frequencies for `background = "model"`.

- test:

  A `character` scalar specifying the type of enrichment test to
  perform. One of `"fisher"` (default) or `"binomial"`. The enrichment
  test is one-sided (enriched in foreground).

- includeRevComp:

  A `logical` scalar. If `TRUE` (default), count k-mer occurrences in
  both `seqs` and their reverse-complement, by concatenating `seqs` and
  their reverse-complemented versions before the counting. This is
  useful if motifs can be expected to occur on any strand (e.g. DNA
  sequences of ChIP-seq peaks). If motifs are only expected on the
  forward strand (e.g. RNA sequences of CLIP-seq peaks),
  `includeRevComp = FALSE` should be used. Note that `bins` will be
  recycled for the reverse complemented sequences, which means that each
  reverse-complemented sequence will be assigned to the same bib as the
  corresponding forward sequence.

- maxFracN:

  A numeric scalar with the maximal fraction of N bases allowed in a
  sequence (defaults to 0.7). Sequences with higher fractions are
  excluded from the analysis.

- maxKmerSize:

  The maximum k-mer size to consider, when adjusting background sequence
  weights for k-mer composition compared to the foreground sequences.
  The default value (3) will correct for mono-, di- and tri-mer
  composition.

- GCbreaks:

  The breaks between GC bins. The default value is based on the
  hard-coded bins used in Homer.

- pseudocount.kmers:

  A `numeric` scalar - will be added to the observed and expected counts
  for each k-mer to avoid zero values.

- pseudocount.log2enr:

  A numerical scalar with the pseudocount to add to foreground and
  background counts when calculating log2 motif enrichments

- p.adjust.method:

  A character scalar selecting the p value adjustment method (used in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)).

- genome:

  A `BSgenome` or `DNAStringSet` object with the genome sequence. Only
  used for `background = "genome"` for extracting background sequences.

- genome.regions:

  An optional
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object defining the intervals in `genome` from which background
  sequences are sampled for `background = "genome"`. If `NULL`,
  background sequences are sampled randomly from `genome`.

- genome.oversample:

  A `numeric` scalar of at least 1.0 defining how many background
  sequences will be sampled per foreground sequence for
  `background = "genome"`. Larger values will take longer but improve
  the sequence composition similarity between foreground and background
  (see `"Details"`).

- BPPARAM:

  An optional
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance determining the parallel back-end to be used during
  evaluation.

- verbose:

  A `logical` scalar. If `TRUE`, report on progress.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with motifs in rows and bins in columns, containing seven assays:

- negLog10P:

  : -log10 P values

- negLog10Padj:

  : -log10 adjusted P values

- pearsonResid:

  : k-mer enrichments as Pearson residuals

- expForegroundWgtWithHits:

  : expected number of foreground sequences with motif hits

- log2enr:

  : k-mer enrichments as log2 ratios

- sumForegroundWgtWithHits:

  : Sum of foreground sequence weights in a bin that have k-mer
  occurrences

- sumBackgroundWgtWithHits:

  : Sum of background sequence weights in a bin that have k-mer
  occurrences

\#' The `rowData` of the object contains annotations (name, PFMs, PWMs
and GC fraction) for the k-mers, while the `colData` slot contains
summary information about the bins.

## Details

This function implements a binned k-mer enrichment analysis. In each
enrichment analysis, the sequences in a specific bin are used as
foreground sequences to test for k-mer enrichments comparing to
background sequences (defined by `background`, see below), similarly as
in done for motifs in
[`calcBinnedMotifEnrR`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedMotifEnrR.md).
Sequences are weighted to correct for GC and shorter k-mer composition
differences between fore- and background sets.

The background sequences are defined according to the value of the
`background` argument:

- otherBins:

  : sequences from all other bins (excluding the current bin)

- allBins:

  : sequences from all bins (including the current bin)

- zeroBin:

  : sequences from the "zero bin", defined by the `maxAbsX` argument of
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md). If
  `bins` does not define a "zero bin", for example because it was
  created by `bin(..., maxAbsX = NULL)`, selecting this background
  definition will abort with an error.

- genome:

  : sequences randomly sampled from the genome (or the intervals defined
  in `genome.regions` if given). For each foreground sequence,
  `genome.oversample` background sequences of the same size are sampled
  (on average). From these, one per foreground sequence is selected
  trying to match the G+C composition. In order to make the sampling
  deterministic, a seed number needs to be provided to the `RNGseed`
  parameter in
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  or
  [`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  when creating the `BiocParallelParam` instance in `BPPARAM`.

- model:

  : a Markov model of the order `MMorder` is estimated from the
  foreground sequences and used to estimate expected k-mer frequencies.
  K-mer enrichments are then calculated comparing observed to these
  expected frequencies. In order to make the process deterministic, a
  seed number needs to be provided to the `RNGseed` parameter in
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  or
  [`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  when creating the `BiocParallelParam` instance in `BPPARAM`.

For each k-mer, the weights of sequences is multiplied with the number
of k-mer occurrences in each sequence and summed, separately for
foreground (`sumForegroundWgtWithHits`) and background
(`sumBackgroundWgtWithHits`) sequences. The function works in ZOOPS
(Zero-Or-One-Per-Sequence) mode, so at most one occurrence per sequence
is counted, which helps reduce the impact of sequence repeats. The total
foreground (`totalWgtForeground`) and background (`totalWgtBackground`)
sum of sequence weights is also calculated. If a k-mer has zero
`sumForegroundWgtWithHits` and `sumBackgroundWgtWithHits`, then any
values (p-values and enrichment) that are calculated using these two
numbers are set to NA.

Two statistical tests for the calculation of enrichment log p-value are
available: `test = "fisher"` (default) to perform Fisher's exact tests,
or `test = "binomial"` to perform binomial tests, using:

- fisher:

  : `fisher.test(x = tab, alternative = "greater")`, where `tab` is the
  contingency table with the summed weights of sequences in foreground
  or background sets (rows), and with or without a occurrences of a
  particular k-mer (columns).

- binomial:

  :
  `pbinom(q = sumForegroundWgtWithHits - 1, size = totalWgtForeground, prob = sumBackgroundWgtWithHits / totalWgtBackground, lower.tail = FALSE, log.p = TRUE)`

## See also

[`getKmerFreq`](https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.md)
used to calculate k-mer enrichments;
[`getSeq,BSgenome-method`](https://rdrr.io/pkg/BSgenome/man/getSeq-methods.html)
which is used to extract sequences from `genomepkg` if `x` is a
`GRanges` object;
[`bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html) that is
used for parallelization;
[`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md) for
binning of regions

## Examples

``` r
seqs <- Biostrings::DNAStringSet(c("GCATGCATGC", "CATGCGCATG"))
bins <- factor(1:2)
calcBinnedKmerEnr(seqs = seqs, bins = bins, kmerLen = 3)
#> class: SummarizedExperiment 
#> dim: 64 2 
#> metadata(5): bins bins.binmode bins.breaks bins.bin0 param
#> assays(7): negLog10P negLog10Padj ... sumForegroundWgtWithHits
#>   sumBackgroundWgtWithHits
#> rownames(64): AAA AAC ... TTG TTT
#> rowData names(5): motif.id motif.name motif.pfm motif.pwm
#>   motif.percentGC
#> colnames(2): 1 2
#> colData names(6): bin.names bin.lower ... totalWgtForeground
#>   totalWgtBackground
```
