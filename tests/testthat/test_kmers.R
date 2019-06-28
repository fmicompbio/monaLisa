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
    expect_is(res1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE),    "data.frame")
    expect_is(res2 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE), "data.frame")
    expect_identical(res1, res2)
    expect_equal(sum(res1$freq.obs), sum(res1$freq.exp), tolerance = 0.001)
    expect_identical(rownames(res1)[which.max(res1$log2enr)], "ACGT")
    expect_true(all(rownames(res1)[res1$FDR < 0.001] %in% c("AACG", "ACGT", "CGTT")))
    ## res1[res1$FDR < 0.001,]

    ## zoops = TRUE
    set.seed(2)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    seqs <- c(DNAStringSet(DNAString(paste(rep("CA", 250), collapse = ""))), seqs[-1])
    expect_is(res3 <- getKmerFreq(seqs, kmerLen = 4, zoops = TRUE),    "data.frame")
    expect_is(res4 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE), "data.frame")
    expect_equal(sum(res3$FDR < 0.001), 0L)
    expect_equal(sum(res4$FDR < 0.001), 2L)
    expect_equal(rownames(res4[res4$FDR < 0.001,]), c("ACAC", "CACA"))
})

test_that("kmerEnrichments works as expected", {
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(SummarizedExperiment)

    gr <- GRanges("chr1", IRanges(sample(2e8, 1000), width = 500))
    seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
    b <- bin(rep(1:5, each = 200), binmode = "equalN", nElements = 200)

    expect_error(kmerEnrichments("error"))
    expect_error(suppressWarnings(kmerEnrichments(gr, b, genomepkg = "not-exising")))
    expect_error(kmerEnrichments(seqs, b[1:20]))
    expect_error(kmerEnrichments(seqs, b, Ncpu = "error"))
    expect_error(kmerEnrichments(seqs, b, Ncpu = 1:2))
    expect_error(kmerEnrichments(seqs, b, Ncpu = -1))

    expect_message(res1 <- kmerEnrichments(as.character(seqs), b, verbose = TRUE))
    res2 <- kmerEnrichments(gr, b, genomepkg = "BSgenome.Hsapiens.UCSC.hg19", verbose = FALSE)
    res3 <- kmerEnrichments(seqs, b, verbose = FALSE)
    res4 <- kmerEnrichments(seqs, as.numeric(b), verbose = FALSE)

    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_is(res3, "SummarizedExperiment")
    expect_is(res4, "SummarizedExperiment")
    expect_identical(assays(res1), assays(res2))
    expect_identical(assays(res1), assays(res3))
    expect_equal(assays(res1), assays(res4), check.attributes = FALSE)
    expect_identical(nrow(rowData(res1)), 1024L)
    expect_identical(nrow(colData(res1)), nlevels(b))
    expect_identical(assayNames(res1), c("p", "FDR", "enr", "log2enr"))
    expect_identical(dim(assay(res1, "log2enr")), c(1024L, 5L))
})
