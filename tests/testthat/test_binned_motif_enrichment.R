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
    bseq2 <- c(b1  = "AGGGGGGGGG", b2  = "AAGGGGGGGG", b3  = "AAAGGGGGGG",
               b4  = "AAAAGGGGGCGC", b5  = "AAAAAGGGGCGCGC", b6  = "AAAAAAGGGG",
               b7  = "AAAAAAAGGG", b8  = "AAAAAAAAGG", b9  = "AAAAAAAAAG",
               b2b = "AAGGGGGGGG", b4b = "AAAAGGGGGG", b6b = "AAAAAAGGGGCGC")
    df <- DataFrame(seqs = DNAStringSet(c(fseq, bseq)),
                    isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq))),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    df2 <- DataFrame(seqs = DNAStringSet(c(fseq, bseq2)),
                     isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq2))),
                     GCfrac = NA_real_,
                     GCbin = NA_integer_,
                     GCwgt = NA_real_,
                     seqWgt = NA_real_)
    attr(df, "err") <- 0
    attr(df2, "err") <- 0
    
    expect_error(.calculateGCweight("error"), "should be a DataFrame")
    expect_error(.calculateGCweight(df, GCbreaks = "error"), "numeric")
    expect_error(.calculateGCweight(df, GCbreaks = 0.2), "length 2 or greater")
    expect_error(.calculateGCweight(df, verbose = "error"), "logical")
    expect_message(.calculateGCweight(df, verbose = TRUE))
    expect_message(.calculateGCweight(df2, verbose = TRUE))

    expect_is(res1 <- .calculateGCweight(df, verbose = FALSE), "DataFrame")
    expect_is(res2 <- .calculateGCweight(df2, verbose = FALSE), "DataFrame")
    expect_is(res3 <- .calculateGCweight(df2, verbose = FALSE, 
                                         normalizeByLength = FALSE), "DataFrame")
    expect_identical(df[-13, 1:2], res1[, 1:2])
    expect_identical(res1$GCfrac, c(9:6,4:1,9:6,4:1,8,6,4) / 10)
    expect_identical(res1[res1$isForeground, "GCwgt"],
                     rep(1.0, sum(res1$isForeground)))
    expect_identical(res1[!res1$isForeground, "GCwgt"],
                     rep(c(1.03125, 0.68750, 1.37500, 1.03125, 0.68750), c(3, 2, 3, 1, 2)))
    expect_identical(res2$GCbin, res3$GCbin)
    expect_identical(res2$GCfrac, res3$GCfrac)
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

    expect_message(res1 <- .calcMotifEnrichment(motifHitMatrix = mhits, df = df, test = "binom", verbose = TRUE))
    expect_is(res1, "data.frame")
    expect_identical(rownames(res1), colnames(mhits))
    expect_identical(round(res1$logP, 3), c(-1.341, -0.038, -0.844, 0))
    expect_message(res2 <- .calcMotifEnrichment(motifHitMatrix = mhits, df = df, test = "fisher", verbose = TRUE))
    expect_is(res2, "data.frame")
    expect_identical(res1[, -2], res2[, -2])
    expect_identical(round(res2$logP, 3), c(-0.99, -0.029, 0, 0))
})


test_that("calcBinnedMotifEnrR() works (synthetic data)", {
    set.seed(1)
    len <- 50
    seqschar <- unlist(lapply(seq.int(150), function(i) {
        paste(sample(c("A","C","G","T"), len, replace = TRUE), collapse = "")
    }))
    b <- factor(rep(c(1, 2, 3), each = 50))
    b <- setZeroBin(b, NA)
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
    gnm <- DNAStringSet(unlist(lapply(1:10,
                                      function(i) paste(sample(c("A","C","G","T"),
                                                               1000 - i * 10,
                                                               replace = TRUE),
                                                        collapse = ""))))
    names(gnm) <- paste0("g", seq_along(gnm))
    
    # argument checks
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm,
                                     min.score = 20), "No motif hits")
    expect_error(calcBinnedMotifEnrR(seqs = Biostrings::DNAStringSet(c("NNN","NNN","NNN")),
                                     bins = factor(1:3), pwmL = pwm), "No sequence passed")
    expect_error(calcBinnedMotifEnrR(seqs = Biostrings::DNAStringSet(c("TAAAAAAT","GCACGTAT","GCCCCCCG")),
                                     bins = factor(1:3), pwmL = pwm, min.score = 6),
                 "No sequences remained after the GC weight")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm,
                                     background = "error"),
                 "should be one of")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm,
                                     background = "zeroBin"),
                 "has to define a zero bin")
    b2 <- b
    attr(b2, "bin0") <- NULL
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b2, pwmL = pwm,
                                     background = "zeroBin"),
                 "has to define a zero bin")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = NULL, pwmL = pwm,
                                     background = "genome", genome = "error"),
                 "'genome' must be")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = NULL, pwmL = pwm,
                                     background = "genome", genome = gnm,
                                     genome.regions = "error"),
                 "'genome.regions' must be")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = NULL, pwmL = pwm,
                                     background = "genome", genome = gnm,
                                     genome.regions = GenomicRanges::GRanges("error", IRanges::IRanges(1, 10))),
                 "seqlevels not contained")
    expect_error(expect_warning(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm,
                                     background = "genome", min.score = 6,
                                     genome = DNAStringSet(c(g1 = paste(rep("G", 300), collapse = "")))),
                                "do not match well"),
                 "No motif hits found in any of the genomic background sequences")
    expect_error(calcBinnedMotifEnrR(seqs = "error"), "DNAStringSet")
    expect_error(calcBinnedMotifEnrR(seqs = as.character(seqs)), "DNAStringSet")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = "error"), "factor")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b[-1]), "must be of equal length")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = "error"), "PWMatrixList")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, maxFracN = "error"), "numeric")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, min.score = 6, maxKmerSize = "error"), "integer")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, pseudocount.log2enr = "error"))
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, pseudocount.pearsonResid = "error"))
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, p.adjust.method = "error"))
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, BPPARAM = "error"), "BiocParallelParam")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, verbose = "error"), "logical")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = TFBSTools::Matrix(pwm[[1]])), "PWMatrixList")

    # results correctness
    expect_message(res1 <- calcBinnedMotifEnr(seqs = seqs,
                                              bins = b,
                                              motifs = pwm,
                                              method = "R",
                                              min.score = 6,
                                              test = "binom",
                                              verbose = TRUE))
    attr(b, "breaks") <- c(1:4 - 0.5)
    expect_message(res2 <- calcBinnedMotifEnr(seqs = seqs,
                                              bins = b,
                                              motifs = pwm,
                                              method = "R",
                                              min.score = 6,
                                              test = "fisher",
                                              verbose = TRUE))
    expect_message(res3 <- calcBinnedMotifEnr(seqs = seqs,
                                              bins = b,
                                              motifs = pwm,
                                              method = "R",
                                              min.score = 6,
                                              test = "fisher",
                                              background = "genome",
                                              genome = gnm,
                                              verbose = TRUE))
    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_is(res3, "SummarizedExperiment")
    expect_identical(dim(res1), c(length(pwm), nlevels(b)))
    expect_identical(dim(res1), dim(res1))
    expect_identical(dimnames(res1), list(names(pwm), levels(b)))
    expect_identical(dimnames(res1), dimnames(res1))
    expect_identical(names(metadata(res1)),
                     c("bins", "bins.binmode", "bins.breaks", "bins.bin0", "param"))
    expect_identical(metadata(res1)$param$test, "binom")
    expect_identical(metadata(res2)$param$test, "fisher")
    expect_identical(metadata(res1)$param[-3], metadata(res2)$param[-3])
    expect_identical(metadata(res2)$bins, b)
    expect_identical(assayNames(res1), c("negLog10P", "negLog10Padj", "pearsonResid", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(assayNames(res2), c("negLog10P", "negLog10Padj", "pearsonResid", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(colnames(rowData(res1)), c("motif.id", "motif.name", "motif.pfm", "motif.pwm", "motif.percentGC"))
    expect_identical(rowData(res1), rowData(res2))
    expect_identical(dim(colData(res1)), c(3L, 6L))
    expect_identical(colData(res1)[, -c(2,3)], colData(res2)[, -c(2,3)])
    expect_identical(colData(res2)[, 2], attr(b, "breaks")[-4])
    expect_identical(colData(res2)[, 3], attr(b, "breaks")[-1])
    expect_equal(-pbinom(q = assay(res1, "sumForegroundWgtWithHits")[, 3] - 1,
                         size = res1$totalWgtForeground[3],
                         prob = assay(res1, "sumBackgroundWgtWithHits")[, 3] /
                         res1$totalWgtBackground[3], lower.tail = FALSE, log.p = TRUE),
                 assay(res1, "negLog10P")[, 3])
    expect_identical(round(assay(res1, "negLog10P"), 3),
                     structure(c(25.893, 0, 0, 34.537, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "negLog10Padj"), 3),
                     structure(c(25.416, 0, 0, 33.759, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "pearsonResid"), 3),
                     structure(c(9.271, -3.449, -3.521, 12.433, -3.853, -4.192),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res1, "log2enr"), 3),
                     structure(c(1.374, -1.219, -1.293, 1.673, -1.299, -1.423),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "negLog10P"), 3),
                     structure(c(17.038, 0, 0, 21.361, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "negLog10Padj"), 3),
                     structure(c(16.561, 0, 0, 20.582, 0, 0),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "pearsonResid"), 3),
                     structure(c(9.271, -3.449, -3.521, 12.433, -3.853, -4.192),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "log2enr"), 3),
                     structure(c(1.374, -1.219, -1.293, 1.673, -1.299, -1.423),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
})

