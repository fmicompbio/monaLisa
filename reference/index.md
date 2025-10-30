# Package index

## monaLisa package overview

- [`monaLisa-package`](https://fmicompbio.github.io/monaLisa/reference/monaLisa-package.md)
  [`monaLisa`](https://fmicompbio.github.io/monaLisa/reference/monaLisa-package.md)
  : monaLisa - MOtif aNAlysis with Lisa.

## Binned motif-enrichment analysis

Bin regions of interest and analyse motifs enrichments in bins.

- [`bin()`](https://fmicompbio.github.io/monaLisa/reference/bin.md) :

  Bin elements of `x`.

- [`getZeroBin()`](https://fmicompbio.github.io/monaLisa/reference/getSetZeroBin.md)
  [`setZeroBin()`](https://fmicompbio.github.io/monaLisa/reference/getSetZeroBin.md)
  : Get and set the zero bin manually

- [`calcBinnedMotifEnrR()`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedMotifEnrR.md)
  :

  Binned Motif Enrichment Analysis with `monaLisa`

- [`calcBinnedMotifEnrHomer()`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedMotifEnrHomer.md)
  : Prepare and run HOMER motif enrichment analysis.

- [`plotBinDiagnostics()`](https://fmicompbio.github.io/monaLisa/reference/plotBinDiagnostics.md)
  : Plot diagnostics of binned sequences

- [`plotBinDensity()`](https://fmicompbio.github.io/monaLisa/reference/plotBinDensity.md)
  : Density plot of binned elements.

- [`plotBinHist()`](https://fmicompbio.github.io/monaLisa/reference/plotBinHist.md)
  : Histogram of binned elements.

- [`plotBinScatter()`](https://fmicompbio.github.io/monaLisa/reference/plotBinScatter.md)
  : Scatter plot (xy-plot) of binned elements.

- [`plotMotifHeatmaps()`](https://fmicompbio.github.io/monaLisa/reference/plotMotifHeatmaps.md)
  : Heatmap of motif enrichments.

## K-mer analysis

Functions for analysing word (k-mer) occurrences in sequences.

- [`getKmerFreq()`](https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.md)
  : Calculate observed and expected k-mer frequencies
- [`calcBinnedKmerEnr()`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedKmerEnr.md)
  : Calculate k-mer enrichment in bins of sequences.

## Identify predictive motifs using stability selection

Functions for selecting motifs that explain observed changes using
Stability Selection.

- [`randLassoStabSel()`](https://fmicompbio.github.io/monaLisa/reference/randLassoStabSel.md)
  : Randomized Lasso Stability Selection
- [`plotSelectionProb()`](https://fmicompbio.github.io/monaLisa/reference/plotSelectionProb.md)
  : Plot selection probabilities of predictors
- [`plotStabilityPaths()`](https://fmicompbio.github.io/monaLisa/reference/plotStabilityPaths.md)
  : Plot Stability Paths

## Miscellaneous

Various utilites and helper functions, typically used in the functions
above.

- [`findMotifHits()`](https://fmicompbio.github.io/monaLisa/reference/findMotifHits-methods.md)
  : Find motif matches in sequences.
- [`annoSeqlogo()`](https://fmicompbio.github.io/monaLisa/reference/annoSeqlogo.md)
  : Sequence logo annotation
- [`dumpJaspar()`](https://fmicompbio.github.io/monaLisa/reference/dumpJaspar.md)
  : Dump Jaspar motifs into a HOMER motif file.
- [`homerToPFMatrixList()`](https://fmicompbio.github.io/monaLisa/reference/homerToPFMatrixList.md)
  : Read a HOMER motif file and create a PFMatrixList
- [`findHomer()`](https://fmicompbio.github.io/monaLisa/reference/findHomer.md)
  : Find HOMER script file.
- [`prepareHomer()`](https://fmicompbio.github.io/monaLisa/reference/prepareHomer.md)
  : Prepare input files for HOMER motif enrichment analysis.
- [`parseHomerOutput()`](https://fmicompbio.github.io/monaLisa/reference/parseHomerOutput.md)
  : Load output from HOMER findMotifsGenome.pl into R
- [`motifSimilarity()`](https://fmicompbio.github.io/monaLisa/reference/motifSimilarity.md)
  : Calculate similarities between pairs of motifs.
- [`motifKmerSimilarity()`](https://fmicompbio.github.io/monaLisa/reference/motifKmerSimilarity.md)
  : Calculate similarities between motifs and k-mers.
- [`seqLogoGrob()`](https://fmicompbio.github.io/monaLisa/reference/seqLogoGrob.md)
  : Create a simple sequence logo grob.
- [`getColsByBin()`](https://fmicompbio.github.io/monaLisa/reference/getColsByBin.md)
  : Get colors by bin.
- [`sampleRandomRegions()`](https://fmicompbio.github.io/monaLisa/reference/sampleRandomRegions.md)
  : Sample random regions of fixed length.
