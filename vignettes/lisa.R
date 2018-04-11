## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
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
#  homerbin <- findHomer(dirs = "/path/to/look/into")
#  
#  # run analysis
#  # (atac.peaks is a GRanges)
#  resL <- runHomer(gr = atac.peaks, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm10",
#                   outdir = "myresults", motifFile = "jaspar2018.motif", scriptFile = homerbin,
#                   regionsize = "given", Ncpu = 4L)
#  

## ---- session------------------------------------------------------------
sessionInfo(package = "lisa")

