# Finding TFs enriched in differentially methylated regions


# load packages
library(GenomicRanges)
library(SummarizedExperiment)
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)


# load regions of interest: low-methylated regions (LMRs)
# - LMRs are provided as a GRanges object for the mm10 genome in the monaLisa
#   package at system.file("extdata", "LMRsESNPmerged.gr.rds", package = "monaLisa")
# - they were generated as described in Stadler, Murr, Burger et al., Nature 2011
#   (https://doi.org/10.1038/nature10716). Briefly:
#   * bisulfite-sequencing data from mouse ES and NP cells were aligned to the
#     mouse mm9 genome using QuasR::qAlign with default parameters
#   * CpG methylation states were called using QuasR::qMeth with default parameters
#   * low-methylated regions were identified using a three-state hidden markov model,
#     but can now be more easily be identified using the MethylSeekR package
#     (https://bioconductor.org/packages/MethylSeekR/)
#   * LMRs were lifted over from mm9 to mm10 coordinates, using
#     rtracklayer::liftOver and the chain file at
#     http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
#   * Average fraction of methylated CpGs in identified LMR regions were
#     calculated for ES and NP cells, and the difference of these fractions
#     (NP - ES) was stored in the metadata column "deltaMeth" of the LMR
#     GRanges object.
regfile <- system.file("extdata", "LMRsESNPmerged.gr.rds", package = "monaLisa")
reg <- readRDS(regfile)
reg


# sub-sample to 10,000 regions
set.seed(1)
reg <- reg[ sample(x = length(reg), size = 10000, replace = FALSE) ]


# make regions of interest equally long
reg <- trim(resize(reg, width = median(width(reg)), fix = "center"))


# extract sequences for regions of interest
regseq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, reg)


# bin regions of interest according to methylation change
# - each bin should contain 800 regions
# - methylation changes below 0.3 are considered as no change
bins <- bin(x = reg$deltaMeth, binmode = "equalN", nElement = 800, minAbsX = 0.3)
table(bins)


# check if bins differ in sequence length or composition
plotBinDiagnostics(seqs = regseq, bins = bins, aspect = "length") # all identical by definition
tapply(width(reg), bins, range)
plotBinDiagnostics(seqs = regseq, bins = bins, aspect = "GCfrac")
plotBinDiagnostics(seqs = regseq, bins = bins, aspect = "dinucfreq")
# regions that lose methyaltion in NP cells (NP-ES < 0, first bin) are GC-poorer


# get vertebrate motifs from Jaspar database
pfms <- getMatrixSet(JASPAR2020, list(tax_group = "vertebrates"))
pwms <- toPWM(pfms)
pwms


# performe binned motif enrichment analysis
se <- calcBinnedMotifEnrR(seqs = regseq,
                          bins = bins,
                          pwmL = pwms,
                          BPPARAM = BiocParallel::MulticoreParam(40),
                          verbose = TRUE)
# saveRDS(se, "inst/extdata/results.binned_motif_enrichment_LMRs.rds")


# perform binned 6-mer enrichment analysis
sk <- calcBinnedKmerEnr(seqs = regseq,
                        bins = bins,
                        kmerLen = 6L,
                        BPPARAM = BiocParallel::MulticoreParam(40),
                        verbose = TRUE)
# saveRDS(sk, "inst/extdata/results.binned_6mer_enrichment_LMRs.rds")
