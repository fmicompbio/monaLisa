<img src="vignettes/monaLisa_logo_v1.png" align="right" alt="monaLisa" width="200px"/>

<br>

# `monaLisa`: MOtif aNAlysis with Lisa

<br>

## Overview

`monaLisa` makes use of her father [Homer](http://homer.ucsd.edu/homer/index.html)
to look for enriched motifs in sets of genomic regions, compared to all other regions.

The motifs come from a collection of transcription factor binding site specificities,
such as `JASPAR2018`. The regions could be for example promoters or accessible regions.
The regions are grouped into bins according to a numerical value assigned to each
bin, such as change of expression or accessibility.

Finally, `monaLisa` can also be used to look for motifs in sequences.

Current contributors include:

- [Michael Stadler](https://github.com/mbstadler)
- [Dania Machlab](https://github.com/machlabd)
- [Lukas Burger](https://github.com/LukasBurger)

## Functionality

A minimal example to run `monaLisa` looks like:

```
# load package
library(monaLisa)

# bin regions
# (atac.peaks.change is a numberical vector)
bins <- bin(x = atac.peaks.change, binmode = "equalN", nElement = 400)

# dump motifs into file
motiffile <- dumpJaspar("jaspar2018.motif", pkg = "JASPAR2018")

# find Homer (findMotifsGenome.pl)
homerfile <- findHomer(dirs = "/path/to/look/into")

# run analysis
# (atac.peaks is a GRanges)
resL <- runHomer(gr = atac.peaks, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm10",
                 outdir = "myresults", motifFile = "jaspar2018.motif", homerfile = homerfile,
                 regionsize = "given", Ncpu = 4L)

```

## Software status

Travis CI (Linux): [![Travis CI build status](https://travis-ci.com/fmicompbio/monaLisa.svg?branch=master)](https://travis-ci.com/fmicompbio/monaLisa)  [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/monaLisa/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/monaLisa)

<!--
## Reference
`EISA` has been described in:  

"Analysis of intronic and exonic reads in RNA-seq data characterizes
transcriptional and post-transcriptional regulation."  
Gaidatzis D., Burger L., Florescu M. and Stadler, M.B.  
*Nat Biotechnol.* **2015**; 33(7):722-9.
[PubMed: 26098447](https://www.ncbi.nlm.nih.gov/pubmed/26098447), [doi: 10.1038/nbt.3269](https://doi.org/10.1038/nbt.3269)

## Download from Bioconductor
[QuasR download page](https://bioconductor.org/packages/QuasR/)

## Software status

| Platforms |  OS  | R CMD check | Coverage | 
|:----------------:|:----------------:|:----------------:|:----------------:|
| Travis CI | Linux | [![Travis CI build status](https://travis-ci.com/fmicompbio/QuasR.svg?branch=master)](https://travis-ci.com/fmicompbio/QuasR) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/QuasR/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/QuasR) |
| Bioc ([_devel_](http://bioconductor.org/packages/devel/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/QuasR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/QuasR) | `NA` |
| Bioc ([_release_](http://bioconductor.org/packages/release/bioc/html/QuasR.html)) | Multiple | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/QuasR.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/QuasR) | `NA` |
-->
