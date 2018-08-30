test_that("findMotifHits() works properly", {
    homerbin <- findHomer("homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.8/bin")

    seqs <- Biostrings::DNAStringSet(c(s1 = "CCCCCAAACCCCC", s2 = "AAAGGGGGAAA"))

    sf <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(x = seqs, filepath = sf)

    M <- matrix(rep(c(.7,.1,.1,.1),3), nrow = 4, dimnames = list(c("A","C","G","T")))
    pwm <- TFBSTools::PWMatrix(ID = "mypwm", name = "mypwm", profileMatrix = M)
    pwmL <- TFBSTools::PWMatrixList(pwm)

    tf <- tempfile(fileext = ".motif")
    lisa:::.dumpPWMsToHomer2File(pwmL, tf)

    # method = "homer2"
    if (!is.na(homerbin)) { # only test if homer2 binary was found
        res1 <- findMotifHits(pwm,  sf,   method = "homer2", homerfile = homerbin) # PWMatrix,character
        res2 <- findMotifHits(pwmL, seqs, method = "homer2", homerfile = homerbin) # PWMatrixList,DNAStringSet
        res3 <- findMotifHits(tf,   sf,   method = "homer2", homerfile = homerbin) # character,character

        expect_true(inherits(res1, "GRanges"))
        expect_equal(length(res1), 3)
        expect_equal(res1$pwmname, c("mypwm","mypwm","mypwm"))
        expect_equal(GenomicRanges::start(res1), c(6, 1, 9))
        expect_identical(res1, res2)
        expect_identical(res1, res3)
    }

    # method = "matchPWM"
    res1 <- findMotifHits(pwm,  sf,   min.score = "90%", method = "matchPWM") # PWMatrix,character
    res2 <- findMotifHits(pwmL, seqs, min.score = "90%", method = "matchPWM") # PWMatrixList,DNAStringSet
    res3 <- findMotifHits(tf,   sf,   min.score = "90%", method = "matchPWM") # character,character

    expect_true(inherits(res1, "GRanges"))
    expect_equal(length(res1), 3)
    expect_equal(as.character(res1$pwmname), rep("mypwm",3))
    expect_equal(GenomicRanges::start(res1), c(6, 1, 9))
    expect_identical(res1, res2)
    expect_identical(res1, res3)

    unlink(c(sf, tf))
})

