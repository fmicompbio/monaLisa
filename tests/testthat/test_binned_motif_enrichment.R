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
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, pseudofreq.pearsonResid = "error"))
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, p.adjust.method = "error"))
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, BPPARAM = "error"), "BiocParallelParam")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = pwm, verbose = "error"), "logical")
    expect_error(calcBinnedMotifEnrR(seqs = seqs, bins = b, pwmL = TFBSTools::Matrix(pwm[[1]])), "PWMatrixList")

    # results correctness
    expect_message(res1 <- calcBinnedMotifEnrR(seqs = seqs,
                                               bins = b,
                                               pwmL = pwm,
                                               min.score = 6,
                                               test = "binom",
                                               verbose = TRUE))
    attr(b, "breaks") <- c(1:4 - 0.5)
    expect_message(res2 <- calcBinnedMotifEnrR(seqs = seqs,
                                               bins = b,
                                               pwmL = pwm,
                                               min.score = 6,
                                               test = "fisher",
                                               verbose = TRUE))
    set.seed(42L)
    expect_message(res3 <- calcBinnedMotifEnrR(seqs = seqs,
                                               bins = b,
                                               pwmL = pwm,
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
    expect_identical(assayNames(res1), c("negLog10P", "negLog10Padj", "pearsonResid", "expForegroundWgtWithHits", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(assayNames(res2), c("negLog10P", "negLog10Padj", "pearsonResid", "expForegroundWgtWithHits", "log2enr", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
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
                     structure(c(7.827, -3.847, -3.902, 8.878, -4.541, -5.114),
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
                     structure(c(7.827, -3.847, -3.902, 8.878, -4.541, -5.114),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
    expect_identical(round(assay(res2, "log2enr"), 3),
                     structure(c(1.374, -1.219, -1.293, 1.673, -1.299, -1.423),
                               dim = c(length(pwm), nlevels(b)),
                               dimnames = list(names(pwm), levels(b))))
})

