test_that("findMotifHits() works properly", {
    homerbin <- findHomer("homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.8/bin")

    if (!is.na(homerbin)) { # only test if homer2 binary was found
        seqs <- Biostrings::DNAStringSet(c(s1="CCCCCAAACCCCC", s2="AAAGGGGGAAA"))

        sf <- tempfile(fileext = ".fa")
        Biostrings::writeXStringSet(x = seqs, filepath = sf)

        M <- matrix(rep(c(.7,.1,.1,.1),3), nrow=4, dimnames=list(c("A","C","G","T")))
        pwm <- TFBSTools::PWMatrix(ID="mypwm", name = "mypwm", profileMatrix=M)
        pwmL <- TFBSTools::PWMatrixList(pwm)

        tf <- tempfile(fileext = ".motif")
        lisa:::.dumpPWMsToHomer2File(pwmL, tf)

        res1 <- findMotifHits(pwm,  sf, homerfile = homerbin) # PWMatrix,character
        res2 <- findMotifHits(pwmL, sf, homerfile = homerbin) # PWMatrixList,character
        res3 <- findMotifHits(tf,   sf, homerfile = homerbin) # character,character

        expect_true(inherits(res1, "GRanges"))
        expect_equal(length(res1), 3)
        expect_equal(res1$pwmname, c("mypwm","mypwm","mypwm"))
        expect_equal(GenomicRanges::start(res1), c(1, 9, 6))
        expect_identical(res1, res2)
        expect_identical(res1, res3)

        unlink(c(sf, tf))
    } else {
        expect_true(file.exists(res))
    }
})

