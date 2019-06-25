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


test_that("compareMotifPair works as expected", {
    res1 <- lisa:::compareMotifPair(m1, m2)
    res2 <- lisa:::compareMotifPair(m1, m3)
    res3 <- lisa:::compareMotifPair(m2, m3)

    expect_identical(res1, list(bestScore = 1.0, bestOffset = 0L, bestDirection = "revcomp"))
    expect_identical(res2, list(bestScore = 0.63146007727884323479, bestOffset = 0L, bestDirection = "revcomp"))
    expect_identical(res3, list(bestScore = 0.63146007727884323479, bestOffset = 7L, bestDirection = "forward"))
})

test_that("motifSimilarity works as expected", {
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

    expect_is(res1 <- motifSimilarity(x = pfmL, y = NULL, method = "R"), "matrix")
    expect_is(res2 <- motifSimilarity(x = pfmL, y = NULL, method = "R", Ncpu = 3L), "matrix")
    expect_is(res3 <- motifSimilarity(x = pfmL, y = pfmL, method = "R"), "matrix")
    expect_is(res4 <- motifSimilarity(x = tmpf, y = NULL, method = "R"), "matrix")
    expect_identical(res1[c(4,7,8)], c(1.0, 0.63146007727884323479, 0.63146007727884323479))
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    expect_identical(res1, res4)

    homerfile <- findHomer("compareMotifs.pl", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin")
    if (!is.na(homerfile)) { # only test at home
        expect_is(res5 <- motifSimilarity(x = tmpf, y = NULL, method = "HOMER", homerfile = homerfile), "matrix")
        expect_equal(res1, res5)
    }

    # clean up
    unlink(x = tmpf)
})

