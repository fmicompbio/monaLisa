context("utils_enrichment")

test_that(".binomEnrichmentTests works", {
    res <- .binomEnrichmentTest(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, verbose = TRUE)
    expect_equal(res[1], 
                 log(binom.test(x = 9, n = 20, p = 5/17, 
                                alternative = "greater")$p.value))
    expect_equal(res[2], 
                 log(binom.test(x = 3, n = 20, p = 8/17, 
                                alternative = "greater")$p.value))
})

test_that(".fisherEnrichmentTest works", {
    res <- .fisherEnrichmentTest(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, verbose = TRUE)
    expect_equal(res[1], 
                 log(fisher.test(x = matrix(c(9, 20 - 9, 5, 17 - 5), nrow = 2),  
                                 alternative = "greater")$p.value))
    expect_equal(res[2], 
                 log(fisher.test(x = matrix(c(3, 20 - 3, 8, 17 - 8), nrow = 2), 
                                 alternative = "greater")$p.value))
})

test_that(".calcPearsonResiduals works", {
    res <- .calcPearsonResiduals(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20)
    expect_equal(res[1], 
                 (9 - (9 + 5)/(20 + 17) * 20)/sqrt((9 + 5)/(20 + 17) * 20 * 
                                                       (1 - 20/(20 + 17)) * 
                                                       (1 - (5 + 9)/(20 + 17))))
    expect_equal(res[2], 
                 (3 - (3 + 8)/(20 + 17) * 20)/sqrt((3 + 8)/(20 + 17) * 20 * 
                                                       (1 - 20/(20 + 17)) * 
                                                       (1 - (8 + 3)/(20 + 17))))
    
})

test_that(".calcExpFg works", {
    res <- .calcExpFg(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20)
    expect_equal(res[1], (9 + 5)/(20 + 17) * 20)
    expect_equal(res[2], (3 + 8)/(20 + 17) * 20)
})

test_that(".calcLog2Enr works", {
    res <- .calcLog2Enr(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, pseudocount = 3)
    expect_equal(res[1], log2((9/20 * 17 + 3)/(5/17 * 17 + 3)))
    expect_equal(res[2], log2((3/20 * 17 + 3)/(8/17 * 17 + 3)))
})

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


test_that(".defineBackground() works", {
    set.seed(1)
    seqs <- DNAStringSet(unlist(lapply(1:90,
                                       function(i) paste(sample(c("A","C","G","T"),
                                                                20, replace = TRUE,
                                                                prob = c(.2,.3,.3,.2)),
                                                         collapse = ""))))
    names(seqs) <- paste0("s", seq_along(seqs))
    b1 <- factor(rep(1:3, each = 30))
    b1 <- setZeroBin(b1, 2)
    b2 <- factor(rep(1:2, each = 45))
    b2 <- setZeroBin(b2, 2)
    gnm <- DNAStringSet(unlist(lapply(1:10,
                                      function(i) paste(sample(c("A","C","G","T"),
                                                               1000 - i * 10,
                                                               replace = TRUE),
                                                        collapse = ""))))
    names(gnm) <- paste0("g", seq_along(gnm))
    gnm2 <- DNAStringSet(unlist(lapply(1:10,
                                       function(i) paste(sample(c("A","C","G","T"),
                                                                1000 - i * 10,
                                                                replace = TRUE,
                                                                prob = c(.4, .1, .1, .4)),
                                                         collapse = ""))))
    names(gnm2) <- paste0("g", seq_along(gnm2))
    
    df1 <- .defineBackground(seqs, b1, "otherBins", 1, NULL, NULL, 2, 0.7)
    df2 <- .defineBackground(seqs, b1, "allBins",   1, NULL, NULL, 2, 0.7)
    df3 <- .defineBackground(seqs, b1, "zeroBin",   1, NULL, NULL, 2, 0.7)
    set.seed(123L)
    df4 <- .defineBackground(seqs, b1, "genome",    1, gnm,  NULL, 2, 0.7)
    
    df5 <- .defineBackground(seqs, b2, "otherBins", 1, NULL, NULL, 2, 0.7)
    df6 <- .defineBackground(seqs, b2, "zeroBin",   1, NULL, NULL, 2, 0.7)
    
    set.seed(42L)
    df7 <- .defineBackground(seqs, b1, "genome",    1, gnm,
                             GenomicRanges::GRanges("g3", IRanges::IRanges(1, 970)),
                             2, 0.7)
    expect_warning(df8 <- .defineBackground(seqs, b1, "genome", 1,
                                            gnm2, NULL, 2, 0.7))
    
    expect_is(df1, "DataFrame")
    expect_is(df2, "DataFrame")
    expect_is(df3, "DataFrame")
    expect_is(df4, "DataFrame")
    expect_is(df5, "DataFrame")
    expect_is(df6, "DataFrame")
    expect_is(df7, "DataFrame")
    expect_is(df8, "DataFrame")
    
    expect_identical(attr(df1, "err"), NA)
    
    expect_identical(dim(df1), c( 90L, 6L))
    expect_identical(dim(df2), c(120L, 6L))
    expect_identical(dim(df3), c( 60L, 6L))
    expect_identical(dim(df4), c( 60L, 6L))
    expect_identical(dim(df5), c( 90L, 6L))
    expect_identical(dim(df6), c( 90L, 6L))
    expect_identical(dim(df7), c( 60L, 6L))
    expect_identical(dim(df8), c( 60L, 6L))
    
    expect_identical(df5, df6)
    
    expect_identical(unlist(lapply(as.character(df7$seqs[!df7$isForeground]),
                                   function(gs) {
                                       grep(gs, as.character(gnm))
                                   }), use.names = FALSE), rep(3L, 30))
    
    expect_true(.checkDfValidity(df1))
    expect_true(.checkDfValidity(df2))
    expect_true(.checkDfValidity(df3))
    expect_true(.checkDfValidity(df4))
    expect_true(.checkDfValidity(df5))
    expect_true(.checkDfValidity(df6))
    expect_true(.checkDfValidity(df7))
    
    expect_identical(sum(df1$isForeground), 30L)
    expect_identical(sum(df2$isForeground), 30L)
    expect_identical(sum(df3$isForeground), 30L)
    expect_identical(sum(df4$isForeground), 30L)
    expect_identical(sum(df5$isForeground), 45L)
    expect_identical(sum(df6$isForeground), 45L)
    expect_identical(sum(df7$isForeground), 30L)
    
    expect_identical(rownames(df1), names(seqs))
    expect_identical(rownames(df2), names(seqs)[c(1:30, 1:90)])
    expect_identical(rownames(df3), names(seqs)[c(1:60)])
    expect_identical(rownames(df4)[1:30], names(seqs)[1:30])
    expect_true(all(grepl("^g[0-9]+$", rownames(df4)[31:60])))
    expect_identical(rownames(df5), names(seqs))
    expect_identical(rownames(df6), names(seqs))
    expect_identical(rownames(df7)[1:30], names(seqs)[1:30])
    expect_true(all(grepl("^g[0-9]+$", rownames(df7)[31:60])))
    
    # for background = "genome", does the sampling make the G+C distribution more similar?
    gnm.tiles <- BSgenome::getSeq(gnm,
                                  unlist(GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(gnm),
                                                                   tilewidth = 20)))
    set.seed(42L)
    df4b <- .defineBackground(gnm.tiles, factor(rep(1, length(gnm.tiles))), "otherBins", 1,
                              NULL, NULL, NULL, 0.7)
    GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)
    gnm.gcbin.tab <- unclass(tabulate(findInterval(df4b$GCfrac, GCbreaks, all.inside = TRUE), 9))
    bg.gcbin.tab  <- unclass(tabulate(findInterval(df4$GCfrac[!df4$isForeground], GCbreaks, all.inside = TRUE), 9))
    fg.gcbin.tab  <- unclass(tabulate(findInterval(df4$GCfrac[ df4$isForeground], GCbreaks, all.inside = TRUE), 9))
    expect_true(cor(fg.gcbin.tab, bg.gcbin.tab) > cor(fg.gcbin.tab, gnm.gcbin.tab))
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
    
    expect_message(res0 <- .iterativeNormForKmers(df, maxKmerSize = 2L, verbose = TRUE))
    expect_is(res0, "DataFrame")
    expect_identical(round(attr(res0, "err"), 6), 0.495758)
    expect_identical(dim(df), dim(res0))
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


test_that(".checkIfSeqsAreEqualLength works as expected", {
    x1 <- Biostrings::DNAStringSet(c("AAA", "CCC", "GGGGGG"))
    x2 <- GenomicRanges::GRanges("chr1",
                                 IRanges::IRanges(start = c(1, 10, 20),
                                                  width = c(3,  3,  6)))
    
    expect_warning(.checkIfSeqsAreEqualLength(x = x1),
                   "Not all elements of x1 have the same length")
    expect_warning(.checkIfSeqsAreEqualLength(x = x2),
                   "Not all elements of x2 have the same length")
    
    expect_true(is.null(.checkIfSeqsAreEqualLength(x = x1[1:2])))
    expect_true(is.null(.checkIfSeqsAreEqualLength(x = x2[1:2])))
})

