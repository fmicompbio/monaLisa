context("k-mers")

# create sequences
set.seed(1)
seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"), size = 500L, replace = TRUE), collapse = ""))
r <- sample(500L - 5L, length(seqs))
substr(seqs, start = r, stop = r + 5L) <- "AACGTT"
seqsDSS <- Biostrings::DNAStringSet(seqs)

test_that("getKmerFreq works as expected", {
    expect_error(getKmerFreq(1L))
    expect_error(getKmerFreq(seqs, kmerLen = "error"))
    expect_error(getKmerFreq(seqs, kmerLen = 1:2))
    expect_error(getKmerFreq(seqs, kmerLen = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = "error"))
    expect_error(getKmerFreq(seqsDSS, MMorder = 1:2))
    expect_error(getKmerFreq(seqsDSS, MMorder = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = -1))
    expect_error(getKmerFreq(seqsDSS, MMorder = 4))

    expect_is(res1 <- getKmerFreq(seqs),    "data.frame")
    expect_is(res2 <- getKmerFreq(seqsDSS), "data.frame")
    expect_identical(res1, res2)
    expect_equal(sum(res1$freq.obs), sum(res1$freq.exp), tolerance = 0.001)
    expect_identical(rownames(res1)[which.max(res1$log2enr)], "ACGT")
    expect_true(all(rownames(res1)[res1$FDR < 0.05] %in% c("AACG", "ACGT", "CGTT")))
})

