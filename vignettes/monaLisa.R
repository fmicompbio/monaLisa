## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  out.width = "80%",
  fig.align = "center"
)

## ---- quick, eval=FALSE-------------------------------------------------------
#  # load package
#  library(monaLisa)
#  
#  # bin regions
#  # (atac.peaks.change is a numerical vector)
#  bins <- bin(x = atac.peaks.change, binmode = "equalN", nElement = 400)
#  
#  # dump motifs into file
#  dumpJaspar("jaspar2018.motif", pkg = "JASPAR2018")
#  
#  # find Homer (findMotifsGenome.pl)
#  homerfile <- findHomer(dirs = "/path/to/look/into")
#  
#  # run analysis
#  # (atac.peaks is a GRanges)
#  se <- runHomer(gr = atac.peaks, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm10",
#                 outdir = "myresults", motifFile = "jaspar2018.motif", homerfile = homerfile,
#                 regionsize = "given", Ncpu = 4L)
#  

## ----loadlib, message=FALSE---------------------------------------------------
library(GenomicRanges)
library(SummarizedExperiment)
library(monaLisa)

## ----loadLMRs-----------------------------------------------------------------
lmrfile <- system.file("extdata", "LMRsESNPmerged.gr.rds", package = "monaLisa")
lmr <- readRDS(lmrfile)
lmr

## ----deltameth----------------------------------------------------------------
hist(lmr$deltaMeth, 100, col = "gray", main = "",
     xlab = "Change of methylation (NP - ES)", ylab = "Number of LMRs")

## ----lmrsel-------------------------------------------------------------------
set.seed(1)
lmrsel <- lmr[ sample(x = length(lmr), size = 10000, replace = FALSE) ]

## ----binlmrs------------------------------------------------------------------
bins <- bin(x = lmrsel$deltaMeth, binmode = "equalN", nElement = 800, minAbsX = 0.3)
table(bins)

## ----plotbins-----------------------------------------------------------------
plotBinDensity(lmrsel$deltaMeth, bins, legend = "topleft")

## ----dumpjaspar---------------------------------------------------------------
motiffile <- tempfile(fileext = ".motif")
dumpJaspar(motiffile, pkg = "JASPAR2018")

## ----homerscript--------------------------------------------------------------
homerfile <- findHomer(dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin/")

## ----runhomer, eval=FALSE-----------------------------------------------------
#  outdir <- tempfile(fileext = ".output")
#  se <- runHomer(gr = lmrsel, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm9",
#                 outdir = outdir, motifFile = motiffile, homerfile = homerfile,
#                 regionsize = "given", Ncpu = 20L)

## ----gethomerresults----------------------------------------------------------
se <- readRDS(system.file("extdata", "se.rds", package = "monaLisa"))

## ----summarizedexperiment-----------------------------------------------------
# summary
se
dim(se) # motifs-by-bins

# motif info
rowData(se)
head(rownames(se))

# bin info
colData(se)
head(colnames(se))

# assays: the motif enrichment results
assayNames(se)
assay(se, "log2enr")[1:5, 1:3]

