context("pfm comparison")

# create artifical motifs, PFMatrix and PFMatrixList
m1 <- matrix(data = c(.1,.3,.3,.3,.01,.01,.01,.97,.5,.0,.0,.5), nrow = 4, dimnames = list(c("A","C","G","T"), NULL))
m2 <- cbind(matrix(.25, nrow = 4, ncol = 2), m1[4:1, ][, 3:1])
rownames(m2) <- rownames(m1)
m3 <- cbind(m1[, c(1,1,1)], matrix(.25, nrow = 4, ncol = 5), m2[, c(5,3,3,1)])

pfm1 <- TFBSTools::PFMatrix(ID = "m1", name = "m1", profileMatrix = 100 * m1)
pfm2 <- TFBSTools::PFMatrix(ID = "m2", name = "m2", profileMatrix = 100 * m2)
pfm3 <- TFBSTools::PFMatrix(ID = "m3", name = "m3", profileMatrix = 100 * m3)

pfmL <- TFBSTools::PFMatrixList(pfm1, pfm2, pfm3)

FourMers <- Biostrings::mkAllStrings(c("A","C","G","T"), 4)


test_that(".compareMotifPair works as expected", {
    res1 <- .compareMotifPair(m1, m2)
    res2 <- .compareMotifPair(m1, m3)
    res3 <- .compareMotifPair(m2, m3)

    expect_identical(res1, list(bestScore = 1.0, bestOffset = 0L, bestDirection = "revcomp"))
    expect_identical(res2, list(bestScore = 0.63146007727884323479, bestOffset = 0L, bestDirection = "revcomp"))
    expect_identical(res3, list(bestScore = 0.63146007727884323479, bestOffset = 7L, bestDirection = "forward"))
})

test_that(".compareMotifKmer works as expected", {
    res1 <- .compareMotifKmer(m1, FourMers)
    res2 <- .compareMotifKmer(m2, FourMers)

    expect_is(res1, "list")
    expect_is(res2, "list")
    expect_length(res1, 2L)
    expect_length(res2, 2L)

    expect_length(b1 <- which(res1$bestScore == max(res1$bestScore)), 42L)
    expect_identical(unique(substr(FourMers[b1], start = 1 + res1$bestOffset[b1], stop = 3 + res1$bestOffset[b1])),
                     c("CTA","CTT","GTA","GTT","TTA","TTT"))

    expect_length(b2 <- which(res2$bestScore == max(res2$bestScore)), 42L)
    o2 <- pmin(5 + res2$bestOffset[b2], 4 - res2$bestOffset[b2], 4)
    expect_identical(unique(substr(FourMers[b2], start = 1, stop = o2)),
                     c("AAA","AAC","AAG",
                       "ATAA","ATAC","ATAG","CAAA","CAAC","CAAG","CTAA","CTAC","CTAG",
                       "GAAA","GAAC","GAAG","GTAA","GTAC","GTAG",
                       "TAA","TAC","TAG",
                       "TTAA","TTAC","TTAG"))
})

test_that("motifSimilarity works as expected", {
    # store motifs in file
    tmpf <- tempfile()
    fh <- file(tmpf, "wt")
    for (i in seq_along(pfmL)) {
        cat(sprintf(">%s\t%s\t5\n", TFBSTools::ID(pfmL[[i]]), TFBSTools::name(pfmL[[i]])),
            file = fh, append = TRUE)
            write.table(file = fh, t(TFBSTools::Matrix(pfmL[[i]])), row.names = FALSE, col.names = FALSE,
                        sep = "\t", quote = FALSE, append = TRUE)
    }
    close(fh)

    expect_error(motifSimilarity(numeric(0)))
    expect_error(motifSimilarity("error"))
    
    bpparams <- BiocParallel::SerialParam()
    if (!identical(.Platform$OS.type, "windows"))
        bpparams <- BiocParallel::MulticoreParam(3)

    expect_is(res1 <- motifSimilarity(x = pfmL, y = NULL, method = "R", 
                                      BPPARAM = bpparams, 
                                      verbose = TRUE), "matrix")
    expect_is(res2 <- motifSimilarity(x = pfmL, y = NULL, method = "R", 
                                      BPPARAM = bpparams, 
                                      verbose = TRUE), "matrix")
    expect_is(res3 <- motifSimilarity(x = pfmL, y = pfmL, method = "R", 
                                      BPPARAM = bpparams, 
                                      verbose = TRUE), "matrix")
    expect_is(res4 <- motifSimilarity(x = pfmL, y = pfmL, method = "R", 
                                      verbose = TRUE), "matrix")
    expect_is(res5 <- motifSimilarity(x = tmpf, y = NULL, method = "R", 
                                      verbose = TRUE), "matrix")
    expect_identical(res1[c(4,7,8)], c(1.0, 0.63146007727884323479, 
                                       0.63146007727884323479))
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    expect_identical(res1, res4)
    expect_identical(res1, res5)
    
    homerfile <- findHomer("compareMotifs.pl", dirs = "/Users/runner/work/monaLisa/monaLisa/homer/bin")
    if (is.na(homerfile)) {
        homerfile <- findHomer("compareMotifs.pl", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.11/bin")
    }
    if (!is.na(homerfile)) { # only test at home
        expect_is(res6 <- motifSimilarity(x = tmpf, y = NULL, method = "HOMER", 
                                          homerfile = homerfile, verbose = TRUE), 
                  "matrix")
        expect_equal(res1, res6)
    }

    # clean up
    unlink(x = tmpf)
})

test_that("motifKmerSimilarity works as expected", {
    # store motifs in file
    tmpf <- tempfile()
    fh <- file(tmpf, "wt")
    for (i in seq_along(pfmL)) {
        cat(sprintf(">%s\t%s\t5\n", TFBSTools::ID(pfmL[[i]]), TFBSTools::name(pfmL[[i]])),
            file = fh, append = TRUE)
        write.table(file = fh, t(TFBSTools::Matrix(pfmL[[i]])), 
                    row.names = FALSE, col.names = FALSE,
                    sep = "\t", quote = FALSE, append = TRUE)
    }
    close(fh)

    bpparams <- BiocParallel::SerialParam()
    if (!identical(.Platform$OS.type, "windows"))
        bpparams <- BiocParallel::MulticoreParam(2)
    
    expect_error(motifKmerSimilarity(1L))
    expect_error(motifKmerSimilarity("error"))
    expect_error(motifKmerSimilarity(pfmL, kmerLen = "error"))
    expect_error(motifKmerSimilarity(pfmL, kmerLen = 2:3))
    expect_error(motifKmerSimilarity(pfmL, kmerLen = 2.5))
    expect_error(motifKmerSimilarity(pfmL, kmerLen = -3))

    expect_is(res1 <- motifKmerSimilarity(x = pfmL, kmerLen = 4, 
                                          verbose = TRUE), "matrix")
    expect_is(res2 <- motifKmerSimilarity(
        x = tmpf, kmerLen = 4L, 
        BPPARAM = bpparams, verbose = TRUE), "matrix")
    expect_is(res3 <- motifKmerSimilarity(x = pfmL, kmers = FourMers, 
                                          verbose = TRUE), "matrix")
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    expect_identical(dim(res1), c(3L, 256L))
    expect_identical(dimnames(res1), list(name(pfmL), FourMers))
    expect_identical(which.max(res1), 2L)
    expect_identical(max(res1), prod(m2[1,-1]))
    
    # test with subset of 4-mers
    FourMersSub <- c("AAAA", "GGCG", "TCGT", "ACTA")
    FourMersSubRevComp <- c("TTTT", "CGCC", "ACGA", "TAGT")
    expect_is(res4 <- motifKmerSimilarity(x = pfmL, kmers = FourMersSub,
                                          verbose = TRUE), "matrix")
    expect_identical(dim(res4), c(3L, 4L))
    expect_identical(res1[, FourMersSub], res4)
    
    # test with reverse complements
    expect_is(res5 <- motifKmerSimilarity(x = pfmL, kmers = FourMersSub,
                                          includeRevComp = TRUE, 
                                          verbose = TRUE), "matrix")
    expect_identical(dim(res5), c(3L, 4L))
    expect_identical(pmax(res1[, FourMersSub], res1[, FourMersSubRevComp]),
                     res5)

    # test with a single 4-mer
    expect_is(res6 <- motifKmerSimilarity(x = pfmL, kmers = "ACGT",
                                          verbose = TRUE), "matrix")
    expect_identical(dim(res6), c(3L, 1L))
    expect_identical(res1[, "ACGT", drop = FALSE], res6)
    
    # clean up
    unlink(x = tmpf)
})

