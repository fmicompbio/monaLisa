#
# Here, we make a synthetic genome to use in the example shown in `calcBinnedMotifEnrHomer`.
#
#

basedir <- "/tungstenfs/groups/gbioinfo/machdani/monaLisa_synthetic_genome_creation_for_example"
setwd(basedir)

library(GenomicRanges)
library(TFBSTools)
library(JASPAR2020)
library(Biostrings)
library(GenomeInfoDb)

## generate random sequences for 3 chromosomes
set.seed(42)
genomedir <- file.path(basedir, "genome")
chrsstr <-
    unlist(lapply(1:3, function(i) paste(sample(x = c("A","C","G","T"),
                                                size = 10000,
                                                replace = TRUE,
                                                prob = c(.3,.2,.2,.3)),
                                         collapse = "")))
names(chrsstr) <- paste0("chr", seq_along(chrsstr))
names(chrsstr)
# [1] "chr1" "chr2" "chr3"

## regions of interest as a GRanges object
gr <- GenomicRanges::tileGenome(seqlengths = nchar(chrsstr),
                                tilewidth = 200, cut.last.tile.in.chrom = TRUE)
gr
# GRanges object with 150 ranges and 0 metadata columns:
#         seqnames     ranges strand
#            <Rle>  <IRanges>  <Rle>
#     [1]     chr1      1-200      *
#     [2]     chr1    201-400      *
#     [3]     chr1    401-600      *
#     [4]     chr1    601-800      *
#     [5]     chr1   801-1000      *
#     ...      ...        ...    ...
#   [146]     chr3  9001-9200      *
#   [147]     chr3  9201-9400      *
#   [148]     chr3  9401-9600      *
#   [149]     chr3  9601-9800      *
#   [150]     chr3 9801-10000      *
#   -------
#   seqinfo: 3 sequences from an unspecified genome

## bins
bins <- factor(GenomicRanges::seqnames(gr))
table(bins)
# bins
# chr1 chr2 chr3 
# 50   50   50 

## motifs (to be planted)
selids <- c("MA0139.1", "MA1102.1", "MA0740.1")
pfm <- TFBSTools::getMatrixSet(JASPAR2020, opts = list(ID = selids))
cons <- unlist(lapply(Matrix(pfm), function(x) paste(rownames(x)[apply(x, 2, which.max)], collapse = "")))
cons
#           MA0139.1              MA1102.1              MA0740.1 
# "TGGCCACCAGGGGGCGCTA"      "CACCAGGGGGCACC"      "GGCCACGCCCCCTT" 

## plant motifs
for (chr1 in names(chrsstr)) {
    i <- which(as.character(GenomeInfoDb::seqnames(gr)) == chr1)
    set.seed(42)
    j <- sample(x = GenomicRanges::start(gr)[i], size = round(length(i) / 3))
    m <- match(chr1, names(chrsstr))
    for (j1 in j)
        substring(chrsstr[chr1], first = j1, last = j1 + nchar(cons[m]) - 1) <- cons[m]
}
chrs <- Biostrings::DNAStringSet(chrsstr)
chrs
# DNAStringSet object of length 3:
#     width seq                                                                                                                             names               
# [1] 10000 TGGCCACCAGGGGGCGCTATCACCATTCTCGCTGACAACGTTACTCCGCGTTTGAGGAATGC...GTGAATGAAACCGCCATTTCGCGAACTTACCATCGATAGTGTCCCAGATATTGGAACGGTGG chr1
# [2] 10000 CACCAGGGGGCACCTCTATAACGGCCATTATAAGTGAGTGCAAAACTCCACATGAGGCTTGC...TAAAGGGAGTAAAAACAAGGTAACCTGGATACTAAACTTATTGTGTCTGTCAACACTATGGC chr2
# [3] 10000 GGCCACGCCCCCTTACCGTACTAGTGTCGTATTCAAGATCGACCTAGATCTTGATATAATCA...CGGCAAGACAAATCAACAGTACGCCAACCTTGGATTAGATCTACTAAATGTGACTATTCGGT chr3

## save created genome
dir.create(genomedir)
genomefile <- file.path(genomedir, "exampleGenome.fa")
Biostrings::writeXStringSet(x = chrs, filepath = genomefile, format = "fasta")

## session
date()
# [1] "Thu Sep 16 16:22:12 2021"
sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /tungstenfs/groups/gbioinfo/Appz/easybuild/software/OpenBLAS/0.3.12-GCC-10.2.0/lib/libopenblas_skylakex-r0.3.12.so
# 
# locale:
#     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#     [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] Biostrings_2.60.2    XVector_0.32.0       JASPAR2020_0.99.10   TFBSTools_1.30.0     GenomicRanges_1.44.0 GenomeInfoDb_1.28.4  IRanges_2.26.0      
# [8] S4Vectors_0.30.0     BiocGenerics_0.38.0 
# 
# loaded via a namespace (and not attached):
#     [1] colorspace_2.0-2            rjson_0.2.20                ellipsis_0.3.2              circlize_0.4.13             GlobalOptions_0.1.2        
# [6] clue_0.3-59                 rstudioapi_0.13             bit64_4.0.5                 AnnotationDbi_1.54.1        fansi_0.5.0                
# [11] codetools_0.2-18            splines_4.1.1               R.methodsS3_1.8.1           doParallel_1.0.16           cachem_1.0.6               
# [16] knitr_1.34                  Rsamtools_2.8.0             Cairo_1.5-12.2              seqLogo_1.58.0              annotate_1.70.0            
# [21] cluster_2.1.2               GO.db_3.13.0                png_0.1-7                   R.oo_1.24.0                 readr_2.0.1                
# [26] compiler_4.1.1              httr_1.4.2                  assertthat_0.2.1            Matrix_1.3-4                fastmap_1.1.0              
# [31] limma_3.48.3                htmltools_0.5.2             tools_4.1.1                 igraph_1.2.6                gtable_0.3.0               
# [36] glue_1.4.2                  TFMPvalue_0.0.8             GenomeInfoDbData_1.2.6      reshape2_1.4.4              dplyr_1.0.7                
# [41] Rcpp_1.0.7                  Biobase_2.52.0              vctrs_0.3.8                 rtracklayer_1.52.1          iterators_1.0.13           
# [46] xfun_0.25                   CNEr_1.28.0                 stringr_1.4.0               lifecycle_1.0.0             restfulr_0.0.13            
# [51] poweRlaw_0.70.6             gtools_3.9.2                XML_3.99-0.7                stringdist_0.9.8            zlibbioc_1.38.0            
# [56] scales_1.1.1                BSgenome_1.60.0             hms_1.1.0                   MatrixGenerics_1.4.3        SummarizedExperiment_1.22.0
# [61] RColorBrewer_1.1-2          ComplexHeatmap_2.8.0        yaml_2.2.1                  memoise_2.0.0               ggplot2_3.3.5              
# [66] stabs_0.6-4                 stringi_1.7.4               RSQLite_2.2.8               BiocIO_1.2.0                foreach_1.5.1              
# [71] caTools_1.18.2              BiocParallel_1.26.2         shape_1.4.6                 rlang_0.4.11                pkgconfig_2.0.3            
# [76] matrixStats_0.60.1          bitops_1.0-7                pracma_2.3.3                evaluate_0.14               lattice_0.20-44            
# [81] purrr_0.3.4                 GenomicAlignments_1.28.0    bit_4.0.4                   tidyselect_1.1.1            plyr_1.8.6                 
# [86] magrittr_2.0.1              R6_2.5.1                    generics_0.1.0              monaLisa_0.1.50             DelayedArray_0.18.0        
# [91] DBI_1.1.1                   pillar_1.6.2                survival_3.2-11             KEGGREST_1.32.0             RCurl_1.98-1.4             
# [96] tibble_3.1.4                crayon_1.4.1                utf8_1.2.2                  tzdb_0.1.2                  rmarkdown_2.10             
# [101] GetoptLong_1.0.5            grid_4.1.1                  blob_1.2.2                  digest_0.6.27               xtable_1.8-4               
# [106] R.utils_2.10.1              munsell_0.5.0               glmnet_4.1-2                DirichletMultinomial_1.34.0
