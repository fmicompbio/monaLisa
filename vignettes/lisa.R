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
#  dumpJaspar("jaspar2018.motif", pkg = "JASPAR2018")
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
hist(lmr$deltaMeth, 100, col = "gray", main = "",
     xlab = "Change of methylation (NP - ES)", ylab = "Number of LMRs")

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
                  cluster = TRUE, maxEnr = 2, maxSig = 10)

## ----findMotifs----------------------------------------------------------
# get sequences of promoters as a DNAStringSet
# (could also be a single DNAString, or the name of a fasta file)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
gr <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 1000, downstream = 500)[c(1,4,5,10)]
library(BSgenome.Hsapiens.UCSC.hg19)
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
seqs

# get motifs as a PWMatrixList
# (could also be a single PWMatrix, or the name of a motif file)
library(JASPAR2018)
library(TFBSTools)
pwms <- getMatrixSet(JASPAR2018, list(matrixtype="PWM", tax_group="vertebrates"))
pwms <- pwms[c("MA0885.1","MA0099.3","MA0033.2","MA0037.3","MA0158.1")]
pwms

# predict hits in sequences
res <- findMotifHits(query = pwms, subject = seqs, min.score = 6.0, method = "matchPWM")
res

# ... or using method = "homer2"
homerfile <- findHomer(homerfile = "homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.8/bin/")
res2 <- findMotifHits(query = pwms, subject = seqs, min.score = 6.0, method = "homer2", homerfile = homerfile)
res2

summary(res %in% res2)

# create hit matrix:
# number of site of each motif per sequence
m <- table(seqnames(res), as.character(res$pwmname))
m

## ---- session------------------------------------------------------------
sessionInfo(package = "lisa")

