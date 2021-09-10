# Overview: 
#   Here, we show how the ATAC-seq dataset in monaLisa was generated
#     starting from the peaks as a GRanges object, and the quantified
#     reads across these peaks, per sample. 
#
# Dataset: 
#   The ATAC-seq samples coming from mice at P0 were downloaded from Encode 
#   and consist of two conditions: liver and lung tissues. Each condition had
#   two replicates, and the corresponding bam files were downloaded.
#   The files have the following IDs: 
#     liver: 
#       - replicate 1: ENCFF146ZCO
#       - replicate 2: ENCFF109LQF
#     lung: 
#       - replicate 1: ENCFF203DOC
#       - replicate 2: ENCFF823PTD
#   
# Processed data:
#   Peak calling was then done per condition, using MACS2 (version 2.2.6). Only 
#   autosomal and distal (at least 1kb away from any TSS) peaks were kept, 
#   followed by merging the peaks from both conditions into one set of peaks. 
#   The reads were quantified across the peaks fro each sample, using QuasR
#
# In this script: 
#   Here, we demonstrate how we produced the example dataset used in monaLisa
#   starting from the distal autosomal peaks, and teh corresponding count table
#   across the samples.
#

## wdir
setwd("/tungstenfs/groups/gbioinfo/machdani/vignette_stabsel_dataSet_for_monaLisa_package/monaLisa_inst_scripts/")

## load peaks and cnt table
peaks <- readRDS("../RData/03_autosomal_distal_peaks.rds")
cnt <- readRDS("../RData/03_autosomal_distal_peaks_cnt.rds")

head(peaks)
# GRanges object with 6 ranges and 0 metadata columns:
#            seqnames          ranges strand
#               <Rle>       <IRanges>  <Rle>
#     peak_1     chr1 3119574-3120615      *
#     peak_2     chr1 3416179-3416430      *
#     peak_3     chr1 3455551-3455813      *
#     peak_4     chr1 4077136-4077336      *
#     peak_5     chr1 4077446-4077701      *
#     peak_6     chr1 4282189-4282805      *
#     -------
#     seqinfo: 19 sequences from an unspecified genome; no seqlengths

head(cnt)
#       width liver_rep1 liver_rep2 lung_rep1 lung_rep2
# peak_1  1042         60         29       237       122
# peak_2   252         26         34         8         9
# peak_3   263         26         24         7         8
# peak_4   201          0          1        22        26
# peak_5   256          3          2        28        41
# peak_6   617          6         14       134        96

all(names(peaks)==rownames(cnt))
# [1] TRUE

## calculate average log2 CPMs per condition
cpm <- sweep(x = cnt[, -1], MARGIN = 2, STATS = colSums(cnt[, -1]), FUN = "/")*1e6
log2_cpm <- log2(cpm + 1)
iByS <- split(colnames(cpm), gsub(".{5}$", "", colnames(cpm)))
avr_cpm <- do.call(cbind, lapply(iByS, function(x){rowMeans(cpm[, x, drop=FALSE])}))
log2_avr_cpm <- log2(avr_cpm + 1)
head(log2_avr_cpm)
#           liver      lung
# peak_1 2.5975953 4.1000557
# peak_2 2.2917217 0.8369844
# peak_3 2.0356378 0.7606035
# peak_4 0.1105451 1.6889044
# peak_5 0.3745536 2.0777416
# peak_6 1.2565505 3.5192534

## calculate log2-fold changes
gr <- peaks
gr$logFC_liver_vs_lung <- log2_avr_cpm[, "liver"] -  log2_avr_cpm[, "lung"]
head(gr)
# GRanges object with 6 ranges and 1 metadata column:
#        seqnames          ranges strand | logFC_liver_vs_lung
#           <Rle>       <IRanges>  <Rle> |           <numeric>
# peak_1     chr1 3119574-3120615      * |            -1.50246
# peak_2     chr1 3416179-3416430      * |             1.45474
# peak_3     chr1 3455551-3455813      * |             1.27503
# peak_4     chr1 4077136-4077336      * |            -1.57836
# peak_5     chr1 4077446-4077701      * |            -1.70319
# peak_6     chr1 4282189-4282805      * |            -2.26270
# -------
#     seqinfo: 19 sequences from an unspecified genome; no seqlengths


## randomly sample 10000 peaks (out of 125933)
n <- 10000
set.seed(123)
keep <- sample(x = 1:length(gr), size = n, replace = FALSE)
gr <- gr[keep]


## save for monaLisa
saveRDS(gr, paste0("atac_liver_vs_lung.rds"))


## Session
date()
# [1] "Thu Sep  9 14:13:57 2021"
sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /tungstenfs/groups/gbioinfo/Appz/easybuild/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_skylakex-r0.3.7.so
# 
# locale:
#     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#     [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#     [1] GenomicRanges_1.42.0 GenomeInfoDb_1.26.7  IRanges_2.24.1      
# [4] S4Vectors_0.28.1     BiocGenerics_0.36.1 
# 
# loaded via a namespace (and not attached):
#     [1] zlibbioc_1.36.0        compiler_4.0.5         tools_4.0.5           
# [4] XVector_0.30.0         GenomeInfoDbData_1.2.4 RCurl_1.98-1.3        
# [7] bitops_1.0-7          

