context("k-mers")

test_that(".cons2matrix works as expected", {
    expect_error(.cons2matrix(x = c("error","error")))
    expect_error(.cons2matrix(x = 1L))
    expect_error(.cons2matrix("ACGT", n = "error"))
    
    expect_is(res1 <- .cons2matrix(x = "ACGT", n = 1L), "matrix")
    expect_is(res2 <- .cons2matrix(x = "ACGT", n = 2L), "matrix")
    expect_identical(dim(res1), c(4L, 4L))
    expect_identical(rownames(res1), c("A","C","G","T"))
    expect_equal(res1, diag(4), check.attributes = FALSE)
    expect_identical(res1 * 2L, res2)
})

test_that("getKmerFreq works as expected", {
    library(Biostrings)

    ## truly random sequences...
    set.seed(1)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    r <- sample(500L - 5L, 100L)
    ## ... with a planted 6-mer
    substr(seqs, start = r, stop = r + 5L) <- "AACGTT"
    seqsDSS <- Biostrings::DNAStringSet(seqs)

    expect_error(getKmerFreq(1L))
    expect_error(getKmerFreq(seqs, kmerLen = "error"))
    expect_error(getKmerFreq(seqs, kmerLen = 1:2))
    expect_error(getKmerFreq(seqs, kmerLen = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = "error"))
    expect_error(getKmerFreq(seqsDSS, MMorder = 1:2))
    expect_error(getKmerFreq(seqsDSS, MMorder = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = -1))
    expect_error(getKmerFreq(seqsDSS, MMorder = 4))
    expect_error(getKmerFreq(seqs, 3, includeRevComp = "error"))

    ## zoops = FALSE
    expect_is(res1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE),    "list")
    expect_is(res2 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE), "list")
    expect_identical(res1, res2)
    expect_equal(sum(res1$freq.obs), sum(res1$freq.exp), tolerance = 0.001)
    expect_identical(names(res1$log2enr)[which.max(res1$log2enr)], "ACGT")
    expect_true(all(names(res1$padj)[res1$padj < 0.001] %in% c("AACG", "ACGT", "CGTT")))

    ## zoops = TRUE
    set.seed(2)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    seqs <- c(DNAStringSet(DNAString(paste(rep("CA", 250), collapse = ""))), seqs[-1])
    expect_is(res3 <- getKmerFreq(seqs, kmerLen = 4, zoops = TRUE),  "list")
    expect_is(res4 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE, includeRevComp = FALSE), "list")
    expect_equal(sum(res3$padj < 0.001), 0L)
    expect_equal(sum(res4$padj < 0.001), 2L)
    expect_equal(names(res4$padj[res4$padj < 0.001]), c("ACAC", "CACA"))

    ## strata
    expect_is(res5 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE, strata = 3), "list")
    expect_equal(sum(res1$freq.obs), sum(res5$freq.obs), tolerance = 0.001)
    expect_equal(sum(res1$freq.exp), sum(res5$freq.exp), tolerance = 0.001)
    expect_true(all(res5$strata %in% c(1L, 2L, 3L)))
    expect_length(res5$freq.strata, 3L)
    expect_equal(sum(res5$CpGoe), 204.030341557169, tolerance = 1e-6)
    expect_identical(names(res5$log2enr)[which.max(res5$log2enr)], "ACGT")
    expect_true(all(names(res5$padj)[res5$padj < 0.001] %in% c("AACG", "ACGT", "CGTT")))
})

test_that(".calcKmerEnrichment works", {
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
    k <- 2L

    expect_error(.calcKmerEnrichment(k = "error", df = df),
                 "'k' must be of type 'numeric'")
    expect_error(.calcKmerEnrichment(k = k, df = "error"),
                 "'df' should be a DataFrame")
    expect_error(.calcKmerEnrichment(k = k, df = df,
                                     test = "error"),
                 "should be one of")
    expect_error(.calcKmerEnrichment(k = k, df = df,
                                     test = "fisher", verbose = "error"),
                 "'verbose' must be of type 'logical'")
    
    expect_message(res1 <- .calcKmerEnrichment(k = k, df = df,
                                               test = "binomial", verbose = TRUE))
    expect_message(res3 <- .calcKmerEnrichment(k = k, df = df,
                                               test = "fisher", verbose = TRUE))
    
    expect_is(res1, "data.frame")
    expect_is(res3, "data.frame")
    
    expect_equal(dim(res1), c(4^k, 6))
    expect_identical(res1$sumForegroundWgtWithHits, res3$sumForegroundWgtWithHits)
    expect_identical(res1$sumBackgroundWgtWithHits, res3$sumBackgroundWgtWithHits)
    expect_equal(res1$logP, c(-0.163497930965412, 0, -0.843717559537887, 0, 0,
                              0, 0, 0, 0, 0, -0.843717559537887, 0, 0, 0, 0, 0))
    expect_equal(res3$logP, c(-0.154150679827258, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0))
})


test_that("calcBinnedKmerEnr works as expected", {
    library(SummarizedExperiment)
    library(BiocParallel)
    if (.Platform$OS.type == "unix") {
        pparams <- MulticoreParam(2L)
    } else {
        pparams <- SerialParam()
    }

    set.seed(1)
    k <- 3L
    seqstr <- unlist(lapply(1:200,
                            function(i) paste(sample(c("A","C","G","T"),
                                                     50, replace = TRUE),
                                              collapse = "")))
    b <- bin(rep(1:2, each = 100), binmode = "equalN", nElements = 100)
    m <- structure(c("AAC", "CCA"), names = levels(b))
    mrc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(m)))
    m2 <- rbind(m, mrc)
    rownames(m2) <- NULL
    for (b1 in levels(b))
        substr(seqstr[b == b1], start = 10, stop = 12) <- m[b1]
    seqs <- DNAStringSet(seqstr)
    gnm <- DNAStringSet(unlist(lapply(1:10,
                                      function(i) paste(sample(c("A","C","G","T"),
                                                               1000 - i * 10,
                                                               replace = TRUE),
                                                        collapse = ""))))
    names(gnm) <- paste0("g", seq_along(gnm))
    

    expect_error(calcBinnedKmerEnr("error"))
    expect_error(calcBinnedKmerEnr(seqs, bins = NULL))
    expect_error(calcBinnedKmerEnr(seqs, bins = b[1:20]))
    expect_error(calcBinnedKmerEnr(seqs, bins = as.numeric(b)))
    expect_error(calcBinnedKmerEnr(seqs, b, kmerLen = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "model", MMorder = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "zeroBin"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, test = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, includeRevComp = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, maxFracN = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, maxKmerSize = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, GCbreaks = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudocount.kmers = -1))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudocount.log2enr = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudofreq.pearsonResid = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, zoops = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, p.adjust.method = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm,
                                   genome.regions = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm,
                                   genome.regions = GRanges("error",
                                                            IRanges(1, 10))))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm, genome.oversample = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm, genome.seed = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, BPPARAM = "error"))
    expect_error(calcBinnedKmerEnr(DNAStringSet(rep("NNNNNNNNNN", 10)), background = "model"))

    expect_message(res1 <- calcBinnedKmerEnr(seqs, b, k, includeRevComp = FALSE, verbose = TRUE))
    RNGkind("L'Ecuyer-CMRG")
    set.seed(42L)
    res2 <- calcBinnedKmerEnr(seqs, b, k, background = "genome", genome = gnm,
                              includeRevComp = FALSE, verbose = FALSE, BPPARAM = pparams)
    RNGkind("default")
    res3 <- calcBinnedKmerEnr(seqs, b, k, background = "model",
                              BPPARAM = pparams)
    res4 <- calcBinnedKmerEnr(seqs, b, k, background = "model",
                              test = "binomial")

    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_is(res3, "SummarizedExperiment")
    expect_is(res4, "SummarizedExperiment")
    expect_identical(assayNames(res1),
                     c("negLog10P", "negLog10Padj", "pearsonResid",
                       "expForegroundWgtWithHits", "log2enr",
                       "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(apply(assay(res1, "negLog10Padj"), 2,
                           function(x) names(x)[x > 3]), m)
    expect_identical(apply(assay(res2, "negLog10Padj"), 2,
                           function(x) names(x)[x > 3]), m)
    expect_identical(apply(assay(res3, "negLog10Padj"), 2,
                           function(x) names(x)[x > 0.5]), m2)
    expect_identical(apply(assay(res4, "negLog10Padj"), 2,
                           function(x) names(x)[x > 1]), m2)
    expect_equal(colSums(assay(res1, "negLog10P")),
                 c(`[1,1.5]` = 33.1936891496806, `(1.5,2]` = 31.5395993919718))
    expect_equal(colSums(assay(res3, "negLog10P")),
                 c(`[1,1.5]` = 27.2409566382244, `(1.5,2]` = 23.5971145106083))
    expect_equal(colSums(assay(res4, "negLog10P")),
                 c(`[1,1.5]` = 37.2309483817477, `(1.5,2]` = 30.1142668663514))
    expect_identical(nrow(res1), as.integer(4^k))
    expect_identical(ncol(res1), nlevels(b))
    expect_identical(assay(res1, "sumForegroundWgtWithHits")[,1],
                     assay(res2, "sumForegroundWgtWithHits")[,1])
    expect_identical(assays(res3)[-c(1, 2)], assays(res4)[-c(1, 2)])
    
    skip_on_os("windows")
    ## Note that all tests after this point (within the test_that block) will 
    ## be skipped on windows.
    expect_equal(colSums(assay(res2, "negLog10P")),
                 c(`[1,1.5]` = 24.5599369371223, `(1.5,2]` = 41.9051486968493))
})

