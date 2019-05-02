context("motif finding")

test_that("findMotifHits() works properly", {
    homerbin <- findHomer("homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin")

    seq1 <- Biostrings::DNAString("CCCCCAAACCCCC")
    seqs <- Biostrings::DNAStringSet(c(s1 = "CCCCCAAACCCCC", s2 = "AAAGGGGGAAA"))

    sf <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(x = seqs, filepath = sf)

    M <- matrix(rep(c(.7,.1,.1,.1),3), nrow = 4, dimnames = list(c("A","C","G","T")))
    pwm <- TFBSTools::PWMatrix(ID = "mypwm", name = "mypwm", profileMatrix = log2(M / 0.25))
    pwmL <- TFBSTools::PWMatrixList(pwm)

    tf <- tempfile(fileext = ".motif")
    lisa:::.dumpPWMsToHomer2File(pwmL, tf)

    tf.pwm <- lisa:::.readPWMsFromHomer2File(tf)

    expect_equal(TFBSTools::Matrix(pwmL[[1]]), TFBSTools::Matrix(tf.pwm[[1]]))

    # method = "matchPWM"
    res1 <- findMotifHits(pwm,  sf,   min.score = "90%", method = "matchPWM") # PWMatrix,character
    res2 <- findMotifHits(pwmL, seqs, min.score = "90%", method = "matchPWM") # PWMatrixList,DNAStringSet
    res3 <- findMotifHits(tf,   sf,   min.score = "90%", method = "matchPWM") # character,character

    expect_true(inherits(res1, "GRanges"))
    expect_equal(length(res1), 3)
    expect_equal(as.character(res1$pwmname), rep("mypwm",3))
    expect_equal(GenomicRanges::start(res1), c(6, 1, 9))
    expect_equal(res1, res2)
    expect_equal(res1, res3)

    # method = "homer2"
    if (!is.na(homerbin)) { # only test if homer2 binary was found
        res4 <- findMotifHits(pwm,  sf,   min.score = "90%",   method = "homer2", homerfile = homerbin) # PWMatrix,character
        res5 <- findMotifHits(pwmL, seqs, min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrixList,DNAStringSet
        res6 <- findMotifHits(tf,   sf,   min.score = "90%",   method = "homer2", homerfile = homerbin) # character,character

        expect_true(inherits(res4, "GRanges"))
        expect_equal(length(res4), 3)
        expect_equal(as.character(res4$pwmname), rep("mypwm",3))
        expect_equal(GenomicRanges::start(res4), c(6, 1, 9))
        expect_equal(res4, res5)
        expect_equal(res4, res6)

        # consistency between "matchPWM" and "homer2"
        expect_equivalent(res1, res4)
    }

    unlink(c(sf, tf))
})

