test_that(".checkDfValidity() works", {
    seqchar <- c(s1 = "GTGCATGCAT", s2 = "ACGTACGTAC")
    df <- DataFrame(seqs = DNAStringSet(seqchar),
                    isForeground = c(TRUE, FALSE),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    
    expect_error(.checkDfValidity("error"), "should be a DataFrame")
    expect_error(.checkDfValidity(df[1:3]), "has to have columns")
    expect_error(.checkDfValidity(DataFrame(seqs = seqchar, df[2:6])), "expected types")
    expect_error(.checkDfValidity(DataFrame(df[1], DataFrame(isForeground = 1:2), df[3:6])), "expected types")
    expect_error(.checkDfValidity(DataFrame(df[1:2], DataFrame(GCfrac = seqchar), df[4:6])), "expected types")
    expect_error(.checkDfValidity(DataFrame(df[1:3], DataFrame(GCbin = seqchar), df[5:6])), "expected types")
    expect_error(.checkDfValidity(DataFrame(df[1:4], DataFrame(GCwgt = seqchar), df[6])), "expected types")
    expect_error(.checkDfValidity(DataFrame(df[1:5], DataFrame(seqWgt = seqchar))), "expected types")
    expect_error(.checkDfValidity(df), "attributes: err")
})


test_that(".filterSeqs() works", {
    seqs <- DNAStringSet(c(s1 = "GTGCATGCATACCA", s2 = "ACGNNNGTAC", s3 = "AAA"))

    expect_error(.filterSeqs("error"), "DNAStringSet")
    expect_error(.filterSeqs(seqs, maxFracN = "error"), "numeric")
    expect_error(.filterSeqs(seqs, minLength = -1), "within \\[0,Inf\\]")
    expect_error(.filterSeqs(seqs, maxLength = 0), "within \\[5,Inf\\]")
    expect_error(.filterSeqs(seqs, verbose = "error"), "logical")
    expect_message(.filterSeqs(seqs, maxFracN = 0.1, verbose = TRUE))
    
    expect_identical(.filterSeqs(seqs), c(TRUE, TRUE, FALSE))
    expect_identical(.filterSeqs(seqs, maxFracN = 0.1, verbose = FALSE), c(TRUE, FALSE, FALSE))
    expect_identical(.filterSeqs(seqs, minLength = 12, verbose = FALSE), c(TRUE, FALSE, FALSE))
    expect_identical(.filterSeqs(seqs, maxLength = 10, verbose = FALSE), c(FALSE, TRUE, FALSE))
})


test_that(".calculateGCweight() works", {
    fseq <- c(f1  = "AGGGGGGGGG", f2  = "AAGGGGGGGG", f3  = "AAAGGGGGGG",
              f4  = "AAAAGGGGGG",                     f6  = "AAAAAAGGGG",
              f7  = "AAAAAAAGGG", f8  = "AAAAAAAAGG", f9  = "AAAAAAAAAG")
    bseq <- c(b1  = "AGGGGGGGGG", b2  = "AAGGGGGGGG", b3  = "AAAGGGGGGG",
              b4  = "AAAAGGGGGG", b5  = "AAAAAGGGGG", b6  = "AAAAAAGGGG",
              b7  = "AAAAAAAGGG", b8  = "AAAAAAAAGG", b9  = "AAAAAAAAAG",
              b2b = "AAGGGGGGGG", b4b = "AAAAGGGGGG", b6b = "AAAAAAGGGG")
    df <- DataFrame(seqs = DNAStringSet(c(fseq, bseq)),
                    isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq))),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    attr(df, "err") <- 0
    
    expect_error(.calculateGCweight("error"), "should be a DataFrame")
    expect_error(.calculateGCweight(df, GCbreaks = "error"), "numeric")
    expect_error(.calculateGCweight(df, GCbreaks = 0.2), "length 2 or greater")
    expect_error(.calculateGCweight(df, verbose = "error"), "logical")
    expect_message(.calculateGCweight(df, verbose = TRUE))

    expect_is(res1 <- .calculateGCweight(df, verbose = FALSE), "DataFrame")
    expect_identical(df[-13, 1:2], res1[, 1:2])
    expect_identical(res1$GCfrac, c(9:6,4:1,9:6,4:1,8,6,4) / 10)
    expect_identical(res1[res1$isForeground, "GCwgt"],
                     rep(1.0, sum(res1$isForeground)))
    expect_identical(res1[!res1$isForeground, "GCwgt"],
                     rep(c(1.03125, 0.68750, 1.37500, 1.03125, 0.68750), c(3, 2, 3, 1, 2)))
})


test_that(".normForKmers() works", {
    set.seed(123)
    len <- 50
    seqschar <- unlist(lapply(seq.int(150), function(i) {
      paste(sample(c("A","C","G","T","N"), len, replace = TRUE, prob = c(.24, .24, .24, .24, .04)), collapse = "")
    }))
    seqs <- DNAStringSet(seqschar)
    kmerfreq <- lapply(1:3, function(k) oligonucleotideFrequency(seqs, width = k))
    gkmers <- lapply(kmerfreq, rowSums)
    kmerfreq <- lapply(1:3, function(k) kmerfreq[[k]] / gkmers[[k]])
    kmerseqrc <- lapply(kmerfreq, function(x) as.character(reverseComplement(DNAStringSet(colnames(x)))))
    seqwgt <- c(rep(1, 75), rnorm(75, mean = 1, sd = 0.1))
    isfg <- rep(c(TRUE, FALSE), each = 75)
    
    res1 <- .normForKmers(kmerfreq, gkmers, kmerseqrc, seqwgt, isfg)

    expect_is(res1, "list")
    expect_length(res1, 2L)
    expect_identical(names(res1), c("seqWgt", "err"))
    expect_identical(round(res1$seqWgt[!isfg], 3),
                     c(0.856, 1.085, 1.06, 0.965, 1.056, 1.011, 1.046, 0.831, 1.204, 
                       0.937, 0.787, 1.101, 0.839, 1.017, 0.882, 0.844, 0.964, 1.04, 
                       0.979, 0.923, 0.967, 1.028, 0.955, 0.92, 1.081, 0.87, 1.064, 
                       0.935, 0.955, 1.141, 0.94, 0.88, 0.873, 1.041, 0.812, 1.023, 
                       0.934, 0.843, 0.986, 1.107, 0.806, 0.876, 1.095, 0.895, 0.866, 
                       0.862, 0.95, 0.987, 0.932, 1.042, 1.049, 1.027, 1.078, 0.98, 
                       0.986, 1.065, 0.994, 1.079, 0.985, 1.023, 1.039, 0.882, 1.029, 
                       0.972, 1.178, 0.931, 1.077, 0.767, 0.923, 0.955, 1.113, 0.949, 
                       0.972, 1.118, 0.989))
})


test_that(".iterativeNormForKmers() works", {
    fseq <- c(f1  = "AGGGGGGGGG", f2  = "AAGGGGGGGG", f3  = "AAAGGGGGGG",
              f4  = "AAAAGGGGGG",                     f6  = "AAAAAAGGGG",
              f7  = "AAAAAAAGGG", f8  = "AAAAAAAAGG", f9  = "AAAAAAAAAG")
    bseq <- c(b1  = "AGGGGGGGGG", b2  = "AAGGGGGGGG", b3  = "AAAGGGGGGG",
              b4  = "AAAAGGGGGG", b5  = "AAAAAGGGGG", b6  = "AAAAAAGGGG",
              b7  = "AAAAAAAGGG", b8  = "AAAAAAAAGG", b9  = "AAAAAAAAAG",
              b2b = "AAGGGGGGGG", b4b = "AAAAGGGGGG", b6b = "AAAAAAGGGG")
    df <- DataFrame(seqs = DNAStringSet(c(fseq, bseq)),
                    isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq))),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    attr(df, "err") <- 0
    df <- .calculateGCweight(df)
    
    expect_error(.iterativeNormForKmers("error"))
    expect_error(.iterativeNormForKmers(df, maxKmerSize = "error"), "integer")
    expect_error(.iterativeNormForKmers(df, minSeqWgt = -1), "within \\(0,Inf\\)")
    expect_error(.iterativeNormForKmers(df, maxIter = "error"), "integer")
    expect_error(.iterativeNormForKmers(df, verbose = "error"), "logical")
    
    expect_message(res1 <- .iterativeNormForKmers(df, verbose = TRUE))
    expect_is(res1, "DataFrame")
    expect_identical(round(attr(res1, "err"), 6), 0.954395)
    attr(df, "err") <- attr(res1, "err")
    expect_identical(dim(res1), dim(df))
    expect_identical(res1[1:5], df[1:5])
    expect_identical(round(res1[!res1$isForeground, "seqWgt"], 3),
                     c(1.251, 0.887, 0.913, 0.675, 0.7, 1.344, 1.347,
                       1.347, 0.887, 0.675, 0.7))
})


test_that(".calcMotifEnrichment works", {
    fseq <- c(f1  = "AGGGGGGGGG", f3  = "AAAGGGGGGG", f4  = "AAAAGGGGGG",
              f6  = "AAAAAAGGGG", f8  = "AAAAAAAAGG")
    bseq <- c(b1  = "AGGGGGGGGG", b2  = "AAGGGGGGGG", b3  = "AAAGGGGGGG",
              b4  = "AAAAGGGGGG", b5  = "AAAAAGGGGG")
    df <- DataFrame(seqs = DNAStringSet(c(fseq, bseq)),
                    isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq))),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    attr(df, "err") <- 0
    df <- .calculateGCweight(df)
    df <- .iterativeNormForKmers(df)
    mhits <- cbind(m1 = c(1, 1, 0, 0, 0, 0, 1),
                   m2 = c(0, 0, 1, 1, 1, 1, 0),
                   m3 = c(1, 1, 1, 1, 1, 1, 1),
                   m4 = c(0, 0, 0, 0, 0, 0, 0))
    
    expect_error(.calcMotifEnrichment("error", df), "has to be a matrix")
    expect_error(.calcMotifEnrichment(mhits, df[1:3, ]), "same number of rows")
    mhits2 <- mhits
    rownames(mhits2) <- as.character(seq.int(nrow(mhits2)))
    expect_error(.calcMotifEnrichment(mhits2, df), "identical rownames")
    expect_error(.calcMotifEnrichment(mhits, df[1:3]), "has to have columns")
    expect_error(.calcMotifEnrichment(mhits, df, test = "error"), "should be one of")
    expect_error(.calcMotifEnrichment(mhits, df, verbose = "error"), "logical")

    # expect_warning(res1 <- .calcMotifEnrichment(motifHitMatrix = mhits, df = df))
    # expect_is(res1, "data.frame")
    expect_is(res1 <- .calcMotifEnrichment(motifHitMatrix = mhits, df = df), "data.frame")
    expect_identical(rownames(res1), colnames(mhits))
    expect_identical(round(res1$logP, 3), c(-1.341, -0.038, -0.844, 0))
    expect_message(res2 <- .calcMotifEnrichment(motifHitMatrix = mhits, df = df, test = "fisher", verbose = TRUE))
    expect_is(res2, "data.frame")
    expect_identical(res1[, -2], res2[, -2])
    expect_identical(round(res2$logP, 3), c(-0.99, -0.029, 0, 0))
})


# test_that("calcBinnedMotifEnrR() works in default mode", {
# 
#     library(GenomicRanges)
#     library(BSgenome.Mmusculus.UCSC.mm10)
#     library(JASPAR2018)
#     library(TFBSTools)
#     library(Biostrings)
#     library(SummarizedExperiment)
#   
#     ## We use pre-selected regions and motifs that we 
#     ## ... know will be enriched in one set of regions vs the other
#     ## ... Below (commented out) is a description of how we chose them:
#     ##
#     ## library(monaLisa)
#     ## library(GenomicRanges)
#     ## library(BSgenome.Mmusculus.UCSC.mm10)
#     ## library(JASPAR2018)
#     ## library(TFBSTools)
#     ## library(Biostrings)
#     ## library(SummarizedExperiment)
#     ## 
#     ## genome <- BSgenome.Mmusculus.UCSC.mm10
#     ##  
#     ## ## Find motif hits on chr10 for all TFs on 100 bp tiles from mm10
#     ## 
#     ## pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates"))
#     ## pwms <- toPWM(pfms) # 579 TFs
#     ##
#     ## seq_info <- seqinfo(genome)
#     ## seq_info <- seq_info["chr10"]
#     ## gr <- tileGenome(seq_info, tilewidth = 100, cut.last.tile.in.chrom = TRUE)
#     ## seqs <- getSeq(genome, gr)
#     ## 
#     ## hits <- findMotifHits(query = pwms, subject = seqs, min.score = 10, method = "matchPWM.concat", BPPARAM = MulticoreParam(10))
#     ## 
#     ## ## Select 2 TFs and the regions (bins) where they have unique hits
#     ## 
#     ## mat <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmname))))
#     ## w <- which(mat > 1)
#     ## mat_zoops <- mat
#     ## mat_zoops[w] <- 1
#     ## 
#     ## keep <- rowSums(mat_zoops) == 1
#     ## mat_zoops_unique <- mat_zoops[keep, ]
#     ## dim(mat_zoops_unique)
#     ## # [1] 19145   574
#     ## 
#     ## mat_zoops_unique[1:6, 1:6]
#     ## #           Ahr::Arnt Alx1 ALX3 Alx4 Ar Arid3a
#     ## # s31072         0    0    0    0  0      0
#     ## # s31095         1    0    0    0  0      0
#     ## # s31096         1    0    0    0  0      0
#     ## # s31097         1    0    0    0  0      0
#     ## # s31098         1    0    0    0  0      0
#     ## # s31100         1    0    0    0  0      0
#     ## 
#     ## head(sort(colSums(mat_zoops_unique), decreasing = TRUE), n=14)
#     ## # ZNF263                RBPJ                SPIB       Nkx2-5(var.2) SMAD2::SMAD3::SMAD4               PRDM1                 MYB 
#     ## # 1055                 534                 463                 394                 346                 337                 305 
#     ## # DMRT3               LIN54               Nr5a2                 VDR               SOX10                Mafb        Pou5f1::Sox2 
#     ## # 286                 285                 263                 263                 257                 251                 245 
#     ##
#     ## ## Choose regions (bins) unique for RBPJ and SPIB respectively 
#     ## RBPJ_seqnames <- rownames(mat_zoops_unique)[as.logical(mat_zoops_unique[, "RBPJ"])]
#     ## SPIB_seqnames <- rownames(mat_zoops_unique)[as.logical(mat_zoops_unique[, "SPIB"])]
#     ## 
#     ## ## Get ranges of bins corresponding to these seqnames (seqname corresponds to index)
#     ## RBPJ_ind <- as.numeric(limma::strsplit2(RBPJ_seqnames, "s")[, 2])
#     ## SPIB_ind <- as.numeric(limma::strsplit2(SPIB_seqnames, "s")[, 2])
#     ## 
#     ## gr_RBPJ <- gr[RBPJ_ind]
#     ## gr_SPIB <- gr[SPIB_ind]
#     ## 
#     ## gr_SPIB_final <- c(gr_SPIB[1:10], gr_RBPJ[1:2])
#     ## gr_RBPJ_final <- c(gr_RBPJ[1:10], gr_SPIB[1:2])
#     ##
#     ## as.data.frame(gr_SPIB_final)
#     ## #      seqnames   start     end width strand
#     ## #   1     chr10 3133401 3133500   100      *
#     ## #   2     chr10 3447301 3447400   100      *
#     ## #   3     chr10 3927801 3927900   100      *
#     ## #   4     chr10 4110201 4110300   100      *
#     ## #   5     chr10 4179501 4179600   100      *
#     ## #   6     chr10 4444901 4445000   100      *
#     ## #   7     chr10 4731101 4731200   100      *
#     ## #   8     chr10 5110901 5111000   100      *
#     ## #   9     chr10 5114401 5114500   100      *
#     ## #   10    chr10 5332501 5332600   100      *
#     ## #   11    chr10 3289801 3289900   100      *
#     ## #   12    chr10 3420201 3420300   100      *
#     ## 
#     ## as.data.frame(gr_RBPJ_final)
#     ## #      seqnames   start     end width strand
#     ## #   1     chr10 3289801 3289900   100      *
#     ## #   2     chr10 3420201 3420300   100      *
#     ## #   3     chr10 3594201 3594300   100      *
#     ## #   4     chr10 3622701 3622800   100      *
#     ## #   5     chr10 3877101 3877200   100      *
#     ## #   6     chr10 3993201 3993300   100      *
#     ## #   7     chr10 4068701 4068800   100      *
#     ## #   8     chr10 4485601 4485700   100      *
#     ## #   9     chr10 4641401 4641500   100      *
#     ## #   10    chr10 4683701 4683800   100      *
#     ## #   11    chr10 3133401 3133500   100      *
#     ## #   12    chr10 3447301 3447400   100      *
#     
#     ## Create dataset
#     gr_SPIB <- GRanges(seqnames = "chr10", 
#                        ranges = IRanges(start = c(3133401, 3447301, 3927801,
#                                                   4110201, 4179501, 4444901,
#                                                   4731101, 5110901, 5114401,
#                                                   5332501, 3289801, 3420201), 
#                                         end = c(3133500, 3447400, 3927900,
#                                                 4110300, 4179600, 4445000,
#                                                 4731200, 5111000, 5114500,
#                                                 5332600, 3289900, 3420300)), 
#                        strand = "*")
#     gr_RBPJ <- GRanges(seqnames = "chr10", 
#                        ranges = IRanges(start = c(3289801, 3420201, 3594201,
#                                                   3622701, 3877101, 3993201,
#                                                   4068701, 4485601, 4641401,
#                                                   4683701, 3133401, 3447301), 
#                                         end = c(3289900, 3420300, 3594300,
#                                                 3622800, 3877200, 3993300,
#                                                 4068800, 4485700, 4641500,
#                                                 4683800, 3133500, 3447400)), 
#                        strand = "*")
#     gr <- c(gr_RBPJ, gr_SPIB)
#     names(gr) <- paste0("peak_", 1:length(gr))
#     seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
#     
#     ## PWMs
#     pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM",
#                                           tax_group = "vertebrates",
#                                           name = c("RBPJ", "SPIB")))
#     pwms <- toPWM(pfms)
#     
#     ## bins
#     b <- bin(x = 1:length(seqs), binmode = "breaks",
#              breaks = c(0, (length(gr_RBPJ)), (length(seqs) + 1)))
#     
#     ## Binned motif enrichment with monaLisa
#     enr_res <- calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwms)
#   
#     ## Tests
#     # ... missing arguments or wrong classes
#     expect_error(calcBinnedMotifEnrR(seqs = "error"), "DNAStringSet")
#     expect_error(calcBinnedMotifEnrR(seqs = as.character(seqs)), "DNAStringSet")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = "error"), "factor")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b[-1]), "must be of equal length")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = "error"), "PWMatrixList")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwms, maxFracN = "error"), "numeric")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwms, maxKmerSize = "error"), "integer")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwms, Ncpu = "error"), "integer")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwms, verbose = "error"), "logical")
#     expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = as.matrix(pwms[[1]])))
# 
#     # ... expected results on dataset
#     expect_is(enr_res, "SummarizedExperiment")
#     expect_identical(rownames(enr_res), c("MA0081.1", "MA1116.1"))
#     expect_identical(unname(round(assay(enr_res, 3), 3)),
#                      matrix(c(-1.930, 2.367, 3.550, -2.238), ncol = 2))
#     expect_identical(unname(round(assay(enr_res, 4), 3)),
#                      matrix(c(-0.436, 0.375, 0.656, -0.729), ncol = 2))
# 
# })

test_that("calcBinnedMotifEnrR() works (synthetic data)", {
    set.seed(1)
    len <- 50
    seqschar <- unlist(lapply(seq.int(150), function(i) {
        paste(sample(c("A","C","G","T"), len, replace = TRUE), collapse = "")
    }))
    b <- factor(rep(c(1, 2, 3), each = 50))
    attr(b, "bin0") <- NA
    m1 <- rbind(A = c(7, 7, 7, 7, 7), # AAAAA
                C = c(1, 1, 1, 1, 1),
                G = c(1, 1, 1, 1, 1),
                T = c(1, 1, 1, 1, 1))
    m2 <- rbind(A = c(7, 1, 1, 1, 7), # ACGTA
                C = c(1, 7, 1, 1, 1),
                G = c(1, 1, 7, 1, 1),
                T = c(1, 1, 1, 7, 1))
    pfm <- TFBSTools::PFMatrixList(m1 = TFBSTools::PFMatrix(ID = "m1", name = "m1", profileMatrix = m1),
                                   m2 = TFBSTools::PFMatrix(ID = "m2", name = "m2", profileMatrix = m2))
    pwm <- TFBSTools::toPWM(pfm)
    # plant m1 in bin 1
    p1 <- sample(len - 2 * ncol(m1), round(sum(b == 1) * 0.8), replace = TRUE)
    substring(seqschar[b == 1], first = p1, last = p1 + ncol(m1) - 1) <- "AAAAA"
    # plant m2 in bin 2
    p2 <- sample(len - 2 * ncol(m2), round(sum(b == 2) * 0.8), replace = TRUE)
    substring(seqschar[b == 2], first = p2, last = p2 + ncol(m2) - 1) <- "ACGTA"
    seqs <- Biostrings::DNAStringSet(seqschar)
    
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm,
                                     min.score = 20), "No motif hits")    
    expect_error(calcBinnedMotifEnrR(seqs = "error"), "DNAStringSet")
    expect_error(calcBinnedMotifEnrR(seqs = as.character(seqs)), "DNAStringSet")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = "error"), "factor")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b[-1]), "must be of equal length")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = "error"), "PWMatrixList")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, maxFracN = "error"), "numeric")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, min.score = 6, maxKmerSize = "error"), "integer")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, BPPARAM = "error"), "BiocParallelParam")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, verbose = "error"), "logical")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = TFBSTools::Matrix(pwm[[1]])), "PWMatrixList")

    expect_message(res1 <- calcBinnedMotifEnr(seqs = seqs,
                                              bins = b,
                                              motifs = pwm,
                                              method = "R",
                                              min.score = 6,
                                              test = "binom",
                                              verbose = TRUE))
    expect_message(res2 <- calcBinnedMotifEnr(seqs = seqs,
                                              bins = b,
                                              motifs = pwm,
                                              method = "R",
                                              min.score = 6,
                                              test = "fisher",
                                              verbose = TRUE))
    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_identical(dim(res1), c(length(pwm), nlevels(b)))
    expect_identical(dim(res1), dim(res1))
    expect_identical(dimnames(res1), list(names(pwm), levels(b)))
    expect_identical(dimnames(res1), dimnames(res1))
    expect_identical(names(metadata(res1)),
                     c("sequences", "bins", "bins.binmode", "bins.breaks", "bins.bin0", 
                       "param.test", "param.maxFracN", "param.maxKmerSize", "param.min.score", 
                       "param.matchMethod", "param.BPPARAM.class", "param.BPARAM.bpnworkers", 
                       "param.verbose"))
    expect_identical(metadata(res1)$param.test, "binom")
    expect_identical(metadata(res2)$param.test, "fisher")
    expect_identical(metadata(res1)[-6], metadata(res2)[-6])
    expect_identical(assayNames(res1), c("p", "FDR", "enr", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(assayNames(res2), c("p", "FDR", "enr", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(colnames(rowData(res1)), c("motif.id", "motif.name", "motif.pfm", "motif.pwm", "motif.percentGC"))
    expect_identical(rowData(res1), rowData(res2))
    expect_identical(dim(colData(res1)), c(3L, 2L))
    expect_identical(colData(res1), colData(res2))
    expect_equal(-pbinom(q = assay(res1, "sumForegroundWgtWithHits")[, 3] - 1,
                         size = se.monaLisa$totalWgtForeground[3],
                         prob = assay(res1, "sumBackgroundWgtWithHits")[, 3] /
                          res1$totalWgtBackground[3], lower.tail = FALSE, log.p = TRUE),
                 assay(res1, "p")[, 3])
    expect_identical(round(assay(res1, "p"), 3),
                     structure(c(25.893, 0, 0, 34.537, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "FDR"), 3),
                     structure(c(25.416, 0, 0, 33.759, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "enr"), 3),
                     structure(c(9.271, -3.449, -3.521, 12.433, -3.853, -4.192),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "log2enr"), 3),
                     structure(c(1.374, -1.219, -1.293, 1.673, -1.299, -1.423),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "p"), 3),
                     structure(c(17.038, 0, 0, 21.361, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "FDR"), 3),
                     structure(c(16.561, 0, 0, 20.582, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "enr"), 3),
                     structure(c(9.271, -3.449, -3.521, 12.433, -3.853, -4.192),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "log2enr"), 3),
                     structure(c(1.374, -1.219, -1.293, 1.673, -1.299, -1.423),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
})

