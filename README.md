<img src="man/figures/monaLisa.png" align="right" alt="monaLisa" width="150"/>

<br>

# `monaLisa`: MOtif aNAlysis with Lisa

<br>

## Overview

`monaLisa` was inspired by her father [Homer](http://homer.ucsd.edu/homer/index.html)
to look for enriched motifs in sets (bins) of genomic regions, compared to all other
regions ("binned motif enrichment analysis").

It uses known motifs representing transcription factor binding preferences,
for example for the `JASPAR2020` Bioconductor package. The regions are for
example promoters or accessible regions, which are grouped into bins according
to a numerical value assigned to each region, such as change of expression
or accessibility. The goal of the analysis is to identify transcription
factors that are associated to that numerical value and thus candidates
to be drivers in the underlying biological process.

In addition to the "binned motif enrichment analysis", `monaLisa` can also be
used to address the above question using stability selection (a form of linear
regression), or to look for motif matches in sequences.

Current contributors include:

- [Michael Stadler](https://github.com/mbstadler)
- [Dania Machlab](https://github.com/machlabd)
- [Lukas Burger](https://github.com/LukasBurger)
- [Charlotte Soneson](https://github.com/csoneson)

## News

- information on the latest changes can be found [here](https://github.com/fmicompbio/monaLisa/blob/master/NEWS.md)
- a preprint is available on [bioRxiv](https://doi.org/10.1101/2021.11.30.470570)
- `monaLisa` is available on [Bioconductor](https://bioconductor.org/packages/monaLisa)
- `monaLisa` is now [published](https://doi.org/10.1093/bioinformatics/btac102)

## Citation

To cite `monaLisa` please use the publication found [here](https://doi.org/10.1093/bioinformatics/btac102) or see `citation("monaLisa")`.

## Installation

`monaLisa` can be installed from Bioconductor via the 
`BiocManager` package:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("monaLisa")
```

## Functionality

Here is a minimal example to run a `monaLisa` analysis:

```
library(monaLisa)

mcparams <- BiocParallel::MulticoreParam(10L)
se <- calcBinnedMotifEnrR(seqs = seqs,   # DNAStringSet (e.g. peak sequences)
                          bins = bins,   # factor that groups 'seqs'
                          pwmL = pwms, # PWMatrixList (know motifs)
                          BPPARAM = mcparams,
                          min.score = 10,
                          verbose = TRUE)
```

The return value `se` is a `SummarizedExperiment` with motifs in rows and bins
in columns, and multiple assays with significance and magnitude of the enrichments.

The inputs for `calcBinnedMotifEnrR` can be easily obtained using other
Bioconductor packages:  
```
# get sequences ('atacPeaks' is a GRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, atacPeaks)

# bin sequences ('atacPeaksChange' is a numerical vector)
bins <- monaLisa::bin(x = atacPeaksChange, binmode = "equalN", nElement = 400)

# obtain known motifs from Jaspar
library(JASPAR2020)
library(TFBSTools)
pwms <- getMatrixSet(JASPAR2020, list(matrixtype = "PWM", tax_group = "vertebrates"))
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

<!-- badges: start -->
Github Actions (multiple OS): [![R build status](https://github.com/fmicompbio/monaLisa/workflows/R-CMD-check/badge.svg)](https://github.com/fmicompbio/monaLisa/actions) [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/monaLisa/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/monaLisa)
<!-- badges: end -->

