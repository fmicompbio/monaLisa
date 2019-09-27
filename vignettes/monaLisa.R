## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  out.width = "80%",
  fig.align = "center"
)

## ---- quick, eval=FALSE----------------------------------------------------
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

## ----loadlib, message=FALSE------------------------------------------------
library(GenomicRanges)
library(SummarizedExperiment)
library(monaLisa)

## ----loadLMRs--------------------------------------------------------------
lmrfile <- system.file("extdata", "LMRsESNPmerged.gr.rds", package = "monaLisa")
lmr <- readRDS(lmrfile)
lmr

## ----deltameth-------------------------------------------------------------
hist(lmr$deltaMeth, 100, col = "gray", main = "",
     xlab = "Change of methylation (NP - ES)", ylab = "Number of LMRs")

## ----lmrsel----------------------------------------------------------------
set.seed(1)
lmrsel <- lmr[ sample(x = length(lmr), size = 10000, replace = FALSE) ]

## ----binlmrs---------------------------------------------------------------
bins <- bin(x = lmrsel$deltaMeth, binmode = "equalN", nElement = 800, minAbsX = 0.3)
table(bins)

## ----plotbins--------------------------------------------------------------
plotBinDensity(lmrsel$deltaMeth, bins, legend = "topleft")

## ----dumpjaspar------------------------------------------------------------
motiffile <- tempfile(fileext = ".motif")
dumpJaspar(motiffile, pkg = "JASPAR2018")

## ----homerscript-----------------------------------------------------------
homerfile <- findHomer(dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin/")

## ----runhomer, eval=FALSE--------------------------------------------------
#  outdir <- tempfile(fileext = ".output")
#  se <- runHomer(gr = lmrsel, b = bins, genomedir = "/work/gbioinfo/DB/genomes/mm9",
#                 outdir = outdir, motifFile = motiffile, homerfile = homerfile,
#                 regionsize = "given", Ncpu = 20L)

## ----gethomerresults-------------------------------------------------------
se <- readRDS(system.file("extdata", "se.rds", package = "monaLisa"))

## ----summarizedexperiment--------------------------------------------------
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

## ----plottfs---------------------------------------------------------------
# select strongly enriched TFs
sel <- apply(assay(se, "log2enr"), 1, function(x) max(abs(x))) > 1.0
sum(sel)
seSel <- se[sel, ]
# shorten names
rownames(seSel) <- sub("\\_.*$","",rownames(seSel))
# plot
plotMotifHeatmaps(x = seSel, which.plots = c("log2enr","FDR"), width = 2.0,
                  cluster = TRUE, maxEnr = 2, maxSig = 10, show_motif_GC = TRUE)

## ----wmclustering, eval=FALSE----------------------------------------------
#  SimMat <- motifSimilarity(rowData(se)$motif.pfm, Ncpu = 20L)

## ----getclusteringresults--------------------------------------------------
SimMat <- readRDS(system.file("extdata", "SimMat.rds", package = "monaLisa"))

## ----checkmatrixorder------------------------------------------------------
all(rownames(SimMat) == rownames(se))

## ----plottfsclustered------------------------------------------------------
# create hclust object, similarity defined by 1 - Pearson correlation
hcl <- hclust(as.dist(1 - SimMat[sel, sel]), method="average")
plotMotifHeatmaps(x = seSel, which.plots = c("log2enr","FDR"), width = 1.2,
                  cluster = hcl, maxEnr = 2, maxSig = 10,
                  show_dendrogram=TRUE, show_seqlogo = TRUE,
                  width.seqlogo = 1.2)

## ----findMotifs------------------------------------------------------------
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
pwms <- getMatrixSet(JASPAR2018, list(matrixtype = "PWM", tax_group = "vertebrates"))
pwms <- pwms[c("MA0885.1","MA0099.3","MA0033.2","MA0037.3","MA0158.1")]
pwms

# predict hits in sequences
res <- findMotifHits(query = pwms, subject = seqs, min.score = 6.0, method = "matchPWM")
res

# ... or using method = "homer2"
homerfile <- findHomer(homerfile = "homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin/")
res2 <- findMotifHits(query = pwms, subject = seqs, min.score = 6.0, method = "homer2", homerfile = homerfile)
res2

summary(res %in% res2)

# create hit matrix:
# number of site of each motif per sequence
m <- table(seqnames(res), as.character(res$pwmname))
m

## ----load_data-------------------------------------------------------------

library(monaLisa)

# Path to extdata 
peaks_path <- system.file("extdata", "Liver_vs_Lung_ATAC_peaks.rds", package = "monaLisa")
response_path <- system.file("extdata", "Liver_vs_Lung_ATAC_logFC.rds", package = "monaLisa")

# Load response vector and peaks GRanges
response <- readRDS(response_path)
peaks <- readRDS(peaks_path)


## ----predictor-------------------------------------------------------------

library(JASPAR2018)
library(TFBSTools)
library("BSgenome.Mmusculus.UCSC.mm10")
library(Biostrings)

# Genome
genome <- BSgenome.Mmusculus.UCSC.mm10

# Get PWMs
pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates")) # for seq-logos
pwms <- toPWM(pfms) # for searching

# Get TFBS on given GRanges
homerfile <- findHomer(homerfile = "homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin/")
hits <- findMotifHits(query = pwms, subject = peaks, min.score = 6.0, method = "homer2",
                      homerfile = homerfile, genome = genome, Ncpu = 2)

# Get predictor matrix
predictor_matrix <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmname))))
predictor_matrix[1:6, 1:6]

# remove peaks that do not have any hits (also from the response vector) 
# and check all peaks in predictor_matrix are in same order as in the 
# peaks variable
w <- names(peaks)[!(names(peaks)%in%rownames(predictor_matrix))]
if(length(w)>0){
  peaks <- peaks[-w]
  response <- response[-w]
}
all(rownames(predictor_matrix)==names(peaks))
all(rownames(predictor_matrix)==names(response))

# Note that in this case all peaks have at least one TFBS for any of the TFs. 
# However, that may not always be the case. If a given set of peaks does 
# not appear in the hits variable one must be careful and remove the missing 
# peaks from the peaks variable before doing the GC content calculations below.

# calculate GC and oeCpG content
peakSeq <- BSgenome::getSeq(genome, peaks)
fMono <- oligonucleotideFrequency(peakSeq, width = 1L, as.prob = TRUE)
fDi <- oligonucleotideFrequency(peakSeq, width = 2L, as.prob = TRUE)
percGC <- fMono[,"G"] + fMono[,"C"]
oeCpG <- (fDi[,"CG"] + 0.01) / (fMono[,"G"] * fMono[,"C"] + 0.01)

# add GC and oeCpG to predictor matrix
predictor_matrix <- cbind(percGC, predictor_matrix)
predictor_matrix <- cbind(oeCpG, predictor_matrix)
predictor_matrix[1:6, 1:6]



## ----run_stability---------------------------------------------------------
library(ComplexHeatmap) # heatmap drawing
library(circlize) # used for color specification

stabs <- randomized_stabsel(x = predictor_matrix, y = response,
                            weakness = 0.8, cutoff = 0.7, PFER = 2, mc.cores = 2)

# plot stability paths ...
plotStabilityPaths(stabs)

# ... and selection probabilities
par(mfrow = c(1,1), mar = c(5,4,4,2) + 3)
plotSelectionProb(stabs)

# plot correlation on TFBS matrix of selected TFs
# ... select predictors for selected TFs and calculate correlation
sel_matrix <- predictor_matrix[, stabs$selected]
sel_cor <- cor(sel_matrix, method = "pearson")

# ... prepare sequence logos
pfmsSel <- pfms[match(rownames(sel_cor), name(pfms))]
maxwidth <- max(sapply(Matrix(pfmsSel), ncol))
seqlogoGrobs <- lapply(pfmsSel, seqLogoGrob, xmax = maxwidth)
hmSeqlogo <- rowAnnotation(logo = anno_seqlogo(seqlogoGrobs, which = "row"),
                           annotation_width = unit(1.5, "inch"),
                           show_legend = FALSE, show_annotation_name = FALSE)

# ... draw heatmap
Heatmap(matrix = sel_cor, name = "Pear. Cor.",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        right_annotation = hmSeqlogo)


## ---- session--------------------------------------------------------------
sessionInfo(package = "monaLisa")

