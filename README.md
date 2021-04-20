<img src="vignettes/monaLisa_logo_v1.png" align="right" alt="monaLisa" width="200px"/>

<br>

# `monaLisa`: MOtif aNAlysis with Lisa

<br>

## Overview

`monaLisa` was inspired by her father [Homer](http://homer.ucsd.edu/homer/index.html)
to look for enriched motifs in sets (bins) of genomic regions, compared to all other
regions ("binned motif enrichment analysis").

It uses known motifs representing transcription factor binding preferences,
for example for the `JASPAR2018` Bioconductor package. The regions are for
example promoters or accessible regions, which are grouped into bins according
to a numerical value assigned to each region, such as change of expression
or accessibility. The goal of the analysis is to identify transcription
factors that are associated to that numerical value and thus candidates
to be drivers in the underlying biological process.

In addition to the "binned motif enrichment analysis", `monaLisa` can also be
used to address the above question using stability selection (a form of linear
regression), or to look for motif matchess in sequences.

Current contributors include:

- [Michael Stadler](https://github.com/mbstadler)
- [Dania Machlab](https://github.com/machlabd)
- [Lukas Burger](https://github.com/LukasBurger)

## Functionality

Here is a minimal example to run a `monaLisa` analysis:

```
library(monaLisa)

mcparams <- BiocParallel::MulticoreParam(10L)
se <- calcBinnedMotifEnr(seqs = seqs,   # DNAStringSet (e.g. peak sequences)
                         bins = bins,   # factor that groups 'seqs'
                         motifs = pwms, # PWMatrixList (know motifs)
                         method = "R",
                         BPPARAM = mcparams,
                         min.score = 10,
                         verbose = TRUE)
```

The return value `se` is a `SummarizedExperiment` with motifs in rows and bins
in columns, and multiple assays with significance and magnitude of the enrichments.

The inputs for `calcBinnedMotifEnr` can be easily obtained using other
Bioconductor packages:  
```
# get sequences ('atacPeaks' is a GRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, atacPeaks)

# bin sequences ('atacPeaksChange' is a numberical vector)
bins <- bin(x = atacPeaksChange, binmode = "equalN", nElement = 400)

# obtain known motifs from Jaspar
library(JASPAR2018)
library(TFBSTools)
pwms <- getMatrixSet(JASPAR2018, list(matrixtype = "PWM", tax_group = "vertebrates"))
```

The results can be conveniently visualized:
```
plotBinDensity(atacPeaksChange, bins, legend = FALSE)
```
<img src="man/figures/monaLisa_binning_small.png" align="center" alt="binning" width="412px"/>

```
plotMotifHeatmaps(se, cluster = TRUE,
                  which.plots = c("enr", "FDR"),
                  show_seqlogo = TRUE)
```
<img src="man/figures/monaLisa_heatmaps_small.png" align="center" alt="heatmaps" width="501px"/>

## Software status

<!-- badges: start -->
Github Actions (multiple OS): [![R build status](https://github.com/fmicompbio/monaLisa/workflows/R-CMD-check/badge.svg)](https://github.com/fmicompbio/monaLisa/actions) [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/monaLisa/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/monaLisa)
<!-- badges: end -->

