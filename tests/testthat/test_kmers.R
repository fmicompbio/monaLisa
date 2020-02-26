context("k-mers")

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

    ## zoops = FALSE
    expect_is(res1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE),    "list")
    expect_is(res2 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE), "list")
    expect_identical(res1, res2)
    expect_equal(sum(res1$freq.obs), sum(res1$freq.exp), tolerance = 0.001)
    expect_identical(names(res1$log2enr)[which.max(res1$log2enr)], "ACGT")
    expect_true(all(names(res1$FDR)[res1$FDR < 0.001] %in% c("AACG", "ACGT", "CGTT")))

    ## zoops = TRUE
    set.seed(2)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    seqs <- c(DNAStringSet(DNAString(paste(rep("CA", 250), collapse = ""))), seqs[-1])
    expect_is(res3 <- getKmerFreq(seqs, kmerLen = 4, zoops = TRUE),  "list")
    expect_is(res4 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE), "list")
    expect_equal(sum(res3$FDR < 0.001), 0L)
    expect_equal(sum(res4$FDR < 0.001), 2L)
    expect_equal(names(res4$FDR[res4$FDR < 0.001]), c("ACAC", "CACA"))

    ## strata
    expect_is(res5 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE, strata = 3), "list")
    expect_equal(sum(res1$freq.obs), sum(res5$freq.obs), tolerance = 0.001)
    expect_equal(sum(res1$freq.exp), sum(res5$freq.exp), tolerance = 0.001)
    expect_true(all(res5$strata %in% c(1L, 2L, 3L)))
    expect_length(res5$freq.strata, 3L)
    expect_equal(sum(res5$CpGoe), 102.0151708, tolerance = 1e-6)
    expect_identical(names(res5$log2enr)[which.max(res5$log2enr)], "ACGT")
    expect_true(all(names(res5$FDR)[res5$FDR < 0.001] %in% c("AACG", "ACGT", "CGTT")))
})

test_that("countKmerPairs works as expected", {
    library(Biostrings)

    seqs <- DNAStringSet(c("AAAAAAN","ATATAT","ACGTAC","N"))

    expect_error(countKmerPairs(x = "error"))
    expect_error(countKmerPairs(x = seqs, k = 0))
    expect_error(countKmerPairs(x = seqs, k = 2, n = 0))

    expect_is(res1 <- countKmerPairs(x = seqs, k = 2, n = 1, zoops = FALSE), "matrix")
    expect_is(res2 <- countKmerPairs(x = seqs, k = 2, n = 1, zoops = TRUE),  "matrix")
    expect_identical(which(res1 > 0, arr.ind = TRUE), which(res2 > 0, arr.ind = TRUE))
    expect_equal(unname(which(res1 > 0, arr.ind = TRUE)),
                 cbind(c(1, 13, 13, 2, 7,  4,  12),
                       c(1,  2,  4, 7, 12, 13, 13)))
    expect_identical(sum(res1), 12)
    expect_identical(sum(res2), 7)
})

test_that("clusterKmers works as expected", {
    library(Biostrings)

    ## truly random sequences...
    set.seed(1)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    r <- sample(500L - 5L, 100L)
    ## ... with a planted 6-mer
    substr(seqs, start = r, stop = r + 5L) <- "AACGTT"
    x1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE)
    x2 <- names(x1$FDR[order(x1$FDR)[1:10]])
    expect_error(clusterKmers(x2, method = "cooccurrence"))
    res1 <- clusterKmers(x1, method = "similarity")
    res2 <- clusterKmers(x2, method = "similarity")
    res3 <- clusterKmers(x2, method = "similarity", allowReverseComplement = TRUE)
    res4 <- clusterKmers(x2, method = "cooccurrence", seqs = DNAStringSet(seqs))
    res5 <- clusterKmers(x2, method = "cooccurrence", seqs = DNAStringSet(seqs), allowReverseComplement = TRUE)

    expect_type(res1, "double")
    expect_length(res1, 10L)
    expect_equal(res1, res2, check.attributes = FALSE)
    expect_identical(as.vector(res1), c(1, 1, 1, 2, 3, 1, 3, 1, 4, 5))
    expect_identical(as.vector(res3), c(1, 1, 1, 2, 3, 1, 3, 1, 3, 4))
    expect_identical(as.vector(res4), c(1, 1, 1, 2, 3, 1, 3, 1, 4, 4))
    expect_identical(as.vector(res5), c(1, 1, 1, 5, 2, 1, 2, 1, 4, 4, 3, 4, 1, 4, 1, 2, 2))
})

test_that("kmerEnrichments works as expected", {
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(SummarizedExperiment)

    set.seed(1)
    gr <- GRanges("chr1", IRanges(sample(2e8, 1000), width = 500))
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
    b <- bin(rep(1:5, each = 200), binmode = "equalN", nElements = 200)

    expect_error(kmerEnrichments("error"))
    expect_error(kmerEnrichments(1L))
    expect_error(suppressWarnings(kmerEnrichments(gr, b, genomepkg = "not-exising")))
    expect_error(kmerEnrichments(seqs, b[1:20]))
    expect_error(kmerEnrichments(seqs, b, background = "error"))
    expect_error(kmerEnrichments(seqs, b, Ncpu = "error"))
    expect_error(kmerEnrichments(seqs, b, Ncpu = 1:2))
    expect_error(kmerEnrichments(seqs, b, Ncpu = -1))

    expect_message(res1 <- kmerEnrichments(as.character(seqs), b, verbose = TRUE))
    res2 <- kmerEnrichments(gr, b, genomepkg = "BSgenome.Hsapiens.UCSC.hg19", verbose = TRUE)
    res3 <- kmerEnrichments(seqs, b, verbose = FALSE)
    res4 <- kmerEnrichments(seqs, as.numeric(b), verbose = TRUE)
    res5 <- kmerEnrichments(seqs, b, background = "model")

    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_is(res3, "SummarizedExperiment")
    expect_is(res4, "SummarizedExperiment")
    expect_is(res5, "SummarizedExperiment")
    expect_identical(assays(res1), assays(res2))
    expect_identical(assays(res1), assays(res3))
    expect_equal(assays(res1), assays(res4), check.attributes = FALSE)
    expect_identical(nrow(rowData(res1)), 1024L)
    expect_identical(nrow(colData(res1)), nlevels(b))
    expect_identical(assayNames(res1), c("p", "FDR", "enr", "log2enr"))
    expect_identical(dim(assay(res1, "log2enr")), c(1024L, 5L))
    expect_equal(diag(cor(assay(res1, "enr"), assay(res5, "enr"))),
                 c(0.4145643727,0.4332920692,0.4948244256,0.3836038418,0.3833744606),
                 check.attributes = FALSE)
})

test_that("convertKmersToMotifs works as expected", {
    library(SummarizedExperiment)
    library(TFBSTools)

    ## truly random sequences...
    nseqs <- 1000
    nbins <- 5
    seqlen <- 100
    set.seed(1)
    seqs <- sapply(seq.int(nseqs), function(i) paste(sample(x = c("A","C","G","T"),
                                                     size = seqlen, replace = TRUE),
                                                     collapse = ""))
    b <- bin(rep(seq.int(nbins), each = round(nseqs / nbins)), binmode = "equalN", nElements = round(nseqs / nbins))
    ## ... with a planted 6-mer
    for (ib in 1:3) {
        i <- which(as.numeric(b) == ib)
        i <- i[seq.int(round(length(i) / ib))]
        r <- sample(seqlen - 5L, length(i), replace = TRUE)
        substr(seqs[i], start = r, stop = r + 5L) <- "AACGTT"
    }
    res1 <- kmerEnrichments(seqs, b, kmerLen = 4, background = "model", verbose = TRUE)
    #o <- order(assay(res1, "log2enr")[, 1], decreasing = TRUE)[1:10]
    #res2 <- plotMotifHeatmaps(res1[o, ], cluster = TRUE, show_dendrogram = TRUE)

    ## motifs...
    pfms <- do.call(PFMatrixList, list(PFMatrix(ID = "m1", name = "m1",
                                                profileMatrix = rbind(A = c(85, 85,  5,  5,  5,  5),
                                                                      C = c( 5,  5, 85,  5,  5,  5),
                                                                      G = c( 5,  5,  5, 85,  5,  5),
                                                                      T = c( 5,  5,  5,  5, 85, 85))),
                                       PFMatrix(ID = "m2", name = "m2",
                                                profileMatrix = rbind(A = c( 5,  5,  5,  5),
                                                                      C = c(85,  5,  5, 85),
                                                                      G = c( 5, 85, 85,  5),
                                                                      T = c( 5,  5,  5,  5)))))
    tf <- tempfile(fileext = ".motif")
    monaLisa:::.dumpPWMsToHomer2File(pwmL = toPWM(pfms), fname = tf)
    pfms <- homerToPFMatrixList(tf) # read back in order to have identical @listData@tags
    a1 <- convertKmersToMotifs(res1, tf, verbose = TRUE)
    a2 <- convertKmersToMotifs(res1, pfms)
    expect_is(a1, "SummarizedExperiment")
    expect_is(a2, "SummarizedExperiment")
    expect_identical(a1, a2)
    expect_equal(rowSums(assay(a1, "enr")), c(m1 = 42.263710609719, m2 = 0.537892231292))
})
