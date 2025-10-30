# Prepare and run HOMER motif enrichment analysis.

Run complete HOMER motif enrichment analysis, consisting of calls to
[`prepareHomer`](https://fmicompbio.github.io/monaLisa/reference/prepareHomer.md),
[`system2`](https://rdrr.io/r/base/system2.html) and
[`parseHomerOutput`](https://fmicompbio.github.io/monaLisa/reference/parseHomerOutput.md).
This function requires `HOMER` to be installed (see
[http://homer.ucsd.edu/homer/index.html](http://homer.ucsd.edu/homer/index.md))
and the path to the tool to be provided (`homerfile` argument).

## Usage

``` r
calcBinnedMotifEnrHomer(
  gr,
  b,
  genomedir,
  outdir,
  motifFile,
  homerfile = findHomer(),
  regionsize = "given",
  pseudocount.log2enr = 8,
  p.adjust.method = "BH",
  Ncpu = 2L,
  verbose = FALSE,
  verbose.Homer = FALSE
)
```

## Arguments

- gr:

  A `GRanges` object (or an object that can be coerced to one) with the
  genomic regions to analyze.

- b:

  A vector of the same length as `gr` that groups its elements into bins
  (typically a factor, such as the one returned by
  [`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md)).

- genomedir:

  Directory containing sequence files in Fasta format (one per
  chromosome).

- outdir:

  A path specifying the folder into which the output files will be
  written.

- motifFile:

  A file with HOMER formatted PWMs to be used in the enrichment
  analysis.

- homerfile:

  Path and file name of the `findMotifsGenome.pl` HOMER script.

- regionsize:

  The peak size to use in HOMER (`"given"` keeps the coordinate region,
  an integer value will keep only that many bases in the region center).

- pseudocount.log2enr:

  A numerical scalar with the pseudocount to add to foreground and
  background counts when calculating log2 motif enrichments

- p.adjust.method:

  A character scalar selecting the p value adjustment method (used in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)).

- Ncpu:

  Number of parallel threads that HOMER can use.

- verbose:

  A logical scalar. If `TRUE`, print progress messages.

- verbose.Homer:

  A logical scalar. If `TRUE`, print the console output when running
  Homer.

## Value

A `SummarizedExperiment` object with motifs in rows and bins in columns,
containing seven assays:

- negLog10P:

  : -log10 P values

- negLog10Padj:

  : -log10 adjusted P values

- pearsonResid:

  : motif enrichments as Pearson residuals

- expForegroundWgtWithHits:

  : expected number of foreground sequences with motif hits

- log2enr:

  : motif enrichments as log2 ratios

- sumForegroundWgtWithHits:

  : Sum of foreground sequence weights in a bin that have motif hits

- sumBackgroundWgtWithHits:

  : Sum of background sequence weights in a bin that have motif hits

The `rowData` of the object contains annotations (name, PFMs, PWMs and
GC fraction) for the motifs, while the `colData` slot contains summary
information about the bins.

## See also

The functions that are wrapped:
[`prepareHomer`](https://fmicompbio.github.io/monaLisa/reference/prepareHomer.md),
[`system2`](https://rdrr.io/r/base/system2.html) and
[`parseHomerOutput`](https://fmicompbio.github.io/monaLisa/reference/parseHomerOutput.md),
[`bin`](https://fmicompbio.github.io/monaLisa/reference/bin.md) for
binning of regions

## Examples

``` r
if (!is.na(findHomer())){

  # genome
  genome <-  system.file("extdata", "exampleGenome.fa", package = "monaLisa")

  # create motif file for Homer
  motiffile <- tempfile()
  motifIDs <- c("MA0139.1", "MA1102.1", "MA0740.1")
  dumpJaspar(filename = motiffile, pkg = "JASPAR2020",
             opts = list(ID = motifIDs))

  # GRanges of regions used in binned motif enrichment analysis
  gr <- GenomicRanges::tileGenome(
      seqlengths = c(chr1 = 10000L, chr2 = 10000L, chr3 = 10000L),
      tilewidth = 200, cut.last.tile.in.chrom = TRUE)

  # create bins (motif enrichment analysis will be per bin)
  bins <- factor(GenomicRanges::seqnames(gr))
  table(bins)

  # run calcBinnedMotifEnrHomer
  outdir <- tempfile()
  se <- calcBinnedMotifEnrHomer(gr = gr, b = bins, genomedir = genome,
      outdir = outdir, motifFile = motiffile)
  list.files(outdir)

  }
```
