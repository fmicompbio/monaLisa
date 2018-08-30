## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  out.width = "80%",
  fig.align = "center"
)

## ---- quick, eval=FALSE--------------------------------------------------
#  # load package
#  library(lisa)
#  
#  # bin regions
#  # (atac.peaks.change is a numberical vector)
#  bins <- bin(x = atac.peaks.change, binmode = "equalN", nElement = 400)
#  
#  # dump motifs into file
#  motiffile <- dumpJaspar("jaspar2018.motif", pkg = "JASPAR2018")
#  
#  # find Homer (findMotifsGenome.pl)
#  homerfile <- findHomer(dirs = "/path/to/look/into")
#  
#  # run analysis
#  # (atac.peaks is a GRanges)
#  resL <- runHomer(gr = atac.peaks, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm10",
#                   outdir = "myresults", motifFile = "jaspar2018.motif", homerfile = homerfile,
#                   regionsize = "given", Ncpu = 4L)
#  

## ----loadlib, message=FALSE----------------------------------------------
library(GenomicRanges)
library(lisa)

## ----loadLMRs------------------------------------------------------------
lmrfile <- system.file("extdata", "LMRsESNPmerged.gr.rds", package = "lisa")
lmr <- readRDS(lmrfile)
lmr

## ----deltameth-----------------------------------------------------------
hist(lmr$deltaMeth, 100, col="gray", main="",
     xlab="Change of methylation (NP - ES)", ylab="Number of LMRs")

## ----lmrsel--------------------------------------------------------------
set.seed(1)
lmrsel <- lmr[ sample(x = length(lmr), size = 10000, replace = FALSE) ]

## ----binlmrs-------------------------------------------------------------
bins <- bin(x = lmrsel$deltaMeth, binmode = "equalN", nElement = 800, minAbsX = 0.3)
table(bins)

## ----plotbins------------------------------------------------------------
plotBinDensity(lmrsel$deltaMeth, bins, legend = "topleft")

## ----dumpjaspar----------------------------------------------------------
motiffile <- tempfile(fileext = ".motif")
dumpJaspar(motiffile, pkg = "JASPAR2018", relScoreCutoff = 0.9)

## ----homerscript---------------------------------------------------------
homerfile <- findHomer(dirs = "/work/gbioinfo/Appz/Homer/Homer-4.8/bin/")

## ----runhomer, eval=FALSE------------------------------------------------
#  outdir <- tempfile(fileext = ".output")
#  resL <- runHomer(gr = lmrsel, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm9",
#                   outdir = outdir, motifFile = motiffile, homerfile = homerfile,
#                   regionsize = "given", Ncpu = 20L)

## ----gethomerresults-----------------------------------------------------
resL <- readRDS(system.file("extdata", "resL.rds", package = "lisa"))

## ----plottfs-------------------------------------------------------------
# select strongly enriched TFs
sel <- apply(resL[["log2enr"]], 1, function(x) max(abs(x))) > 1.0
sum(sel)
resLsel <- lapply(resL, function(x) x[sel,])
# shorten names
resLsel <- lapply(resLsel, function(x) { rownames(x) <- sub("\\|.*$","",rownames(x)); x })
# plot
plotMotifHeatmaps(x = resLsel, b = bins, which.plots = c("log2enr","FDR"), width = 2.0,
                  cluster=TRUE, maxEnr = 2, maxSig = 10)

## ---- session------------------------------------------------------------
sessionInfo(package = "lisa")

