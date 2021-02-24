context("motif finding")

test_that("findMotifHits() works properly", {
    homerbin <- findHomer("homer2", dirs = "/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin")

    seq1 <- Biostrings::DNAString("CCCCCAAACCCCC")
    seqs <- Biostrings::DNAStringSet(c(seq1 = "CCCCCAAACCCCC", seq2 = "AAAGGGGGAAA"))
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(15927933, 3261215), end = c(15927945, 3261225), names = c("seq1", "seq2")), strand = "+")

    sf <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(x = seqs, filepath = sf)

    M <- matrix(rep(c(.7,.1,.1,.1),3), nrow = 4, dimnames = list(c("A","C","G","T")))
    pwm <- TFBSTools::PWMatrix(ID = "mypwm", name = "mypwm", profileMatrix = log2(M / 0.25))
    pwmL <- TFBSTools::PWMatrixList(pwm)
    names(pwmL) <- "mypwm"

    tf <- tempfile(fileext = ".motif")
    monaLisa:::.dumpPWMsToHomer2File(pwmL, tf)

    tf.pwm <- monaLisa:::.readPWMsFromHomer2File(tf)

    expect_equal(TFBSTools::Matrix(pwmL[[1]]), TFBSTools::Matrix(tf.pwm[[1]]))

    # method = "matchPWM"
    ############################################################################
    expect_error(findMotifHits(pwmL, sf,   min.score = TRUE, method = "matchPWM"))
    expect_error(findMotifHits(pwmL, seqs, min.score = TRUE, method = "matchPWM"))

    res1  <- findMotifHits(tf,   sf,   min.score = "90%", method = "matchPWM") # character,character
    res2  <- findMotifHits(tf,   seq1, min.score = "90%", method = "matchPWM") # character,DNAString
    res3  <- findMotifHits(tf,   seqs, min.score = "90%", method = "matchPWM") # character,DNAStringSet
    res4  <- findMotifHits(pwm,  sf,   min.score = "90%", method = "matchPWM") # PWMatrix,character
    res5  <- findMotifHits(pwm,  seq1, min.score = "90%", method = "matchPWM") # PWMatrix,DNAString
    res6  <- findMotifHits(pwm,  seqs, min.score = "90%", method = "matchPWM") # PWMatrix,DNAStringSet
    res7  <- findMotifHits(pwmL, sf,   min.score = "90%", method = "matchPWM") # PWMatrixList,character
    expect_warning(res7b <- findMotifHits(pwmL, sf,   min.score = 4.456, method = "matchPWM")) # PWMatrixList,character
    res8  <- findMotifHits(pwmL, seq1, min.score = "90%", method = "matchPWM") # PWMatrixList,DNAString
    res9  <- findMotifHits(pwmL, seqs, min.score = "90%", method = "matchPWM") # PWMatrixList,DNAStringSet
    expect_warning(res9b <- findMotifHits(pwmL, seqs, min.score = 4.456, method = "matchPWM")) # PWMatrixList,DNAStringSet
    res10 <- findMotifHits(pwm,  gr,   min.score = "90%", method = "matchPWM", genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrix,GRanges
    res11 <- findMotifHits(pwmL, gr,   min.score = "90%", method = "matchPWM", genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrixList,GRanges

    # correctness
    expect_true(inherits(res1, "GRanges"))
    expect_length(res1, 3L)
    expect_equal(as.character(res1$pwmname), rep("mypwm",3))
    expect_equal(GenomicRanges::start(res1), c(6, 1, 9))

    expect_true(inherits(res2, "GRanges"))
    expect_length(res2, 1L)
    expect_equal(as.character(res2$pwmname), "mypwm")
    expect_equal(GenomicRanges::start(res2), 6L)

    # consistency between methods
    expect_equal(res1, res3)
    expect_equal(res1, res4)
    expect_equal(res2, res5)
    expect_equal(res1, res6)
    expect_equal(res1, res7)
    expect_equal(res1, res7b)
    expect_equal(res2, res8)
    expect_equal(res1, res9)
    expect_equal(res1, res9b)
    expect_equal(res1, res10)
    expect_equal(res1, res11)

    # method = "matchPWM.concat"
    ############################################################################
    expect_error(findMotifHits(pwmL, sf,   min.score = TRUE, method = "matchPWM.concat"))
    expect_error(findMotifHits(pwmL, seqs, min.score = TRUE, method = "matchPWM.concat"))
    
    res1c  <- findMotifHits(tf,   sf,   min.score = "90%", method = "matchPWM.concat") # character,character
    res2c  <- findMotifHits(tf,   seq1, min.score = "90%", method = "matchPWM.concat") # character,DNAString
    res3c  <- findMotifHits(tf,   seqs, min.score = "90%", method = "matchPWM.concat") # character,DNAStringSet
    res4c  <- findMotifHits(pwm,  sf,   min.score = "90%", method = "matchPWM.concat") # PWMatrix,character
    res5c  <- findMotifHits(pwm,  seq1, min.score = "90%", method = "matchPWM.concat") # PWMatrix,DNAString
    res6c  <- findMotifHits(pwm,  seqs, min.score = "90%", method = "matchPWM.concat") # PWMatrix,DNAStringSet
    res7c  <- findMotifHits(pwmL, sf,   min.score = "90%", method = "matchPWM.concat") # PWMatrixList,character
    res7cb <- findMotifHits(pwmL, sf,   min.score = 4.456, method = "matchPWM.concat") # PWMatrixList,character
    res8c  <- findMotifHits(pwmL, seq1, min.score = "90%", method = "matchPWM.concat") # PWMatrixList,DNAString
    res9c  <- findMotifHits(pwmL, seqs, min.score = "90%", method = "matchPWM.concat") # PWMatrixList,DNAStringSet
    res9cb <- findMotifHits(pwmL, seqs, min.score = 4.456, method = "matchPWM.concat") # PWMatrixList,DNAStringSet
    res10c <- findMotifHits(pwm,  gr,   min.score = "90%", method = "matchPWM.concat", genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrix,GRanges
    res11c <- findMotifHits(pwmL, gr,   min.score = "90%", method = "matchPWM.concat", genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrixList,GRanges
    
    # correctness
    expect_true(inherits(res1c, "GRanges"))
    expect_length(res1c, 3L)
    expect_equal(as.character(res1c$pwmname), rep("mypwm",3))
    expect_equal(GenomicRanges::start(res1c), c(6, 1, 9))
    
    expect_true(inherits(res2c, "GRanges"))
    expect_length(res2c, 1L)
    expect_equal(as.character(res2c$pwmname), "mypwm")
    expect_equal(GenomicRanges::start(res2c), 6L)
    
    # consistency between "matchPWM" and "matchPWM.concat"
    expect_equivalent(res1, res1c)
    expect_equivalent(res2, res2c)
    
    # consistency between methods
    expect_equal(res1c, res3c)
    expect_equal(res1c, res4c)
    expect_equal(res2c, res5c)
    expect_equal(res1c, res6c)
    expect_equal(res1c, res7c)
    expect_equal(res1c, res7cb)
    expect_equal(res2c, res8c)
    expect_equal(res1c, res9c)
    expect_equal(res1c, res9cb)
    expect_equal(res1c, res10c)
    expect_equal(res1c, res11c)
    
    # method = "homer2"
    if (!is.na(homerbin)) { # only test if homer2 binary was found
        expect_error(findMotifHits(tf, DNAString(paste(rep("A", 1e6), collapse = "")), method = "homer2", homerfile = homerbin))
        expect_error(findMotifHits("not-existing", sf, method = "homer2", homerfile = homerbin))
        expect_error(findMotifHits(tf, "not-existing", method = "homer2", homerfile = homerbin))
        expect_error(findMotifHits(tf, sf, method = "error"))
        expect_error(findMotifHits(pwm,  gr,   min.score = "90%", method = "homer2", homerfile = homerbin))
        expect_error(findMotifHits(pwm,  gr,   min.score = "90%", method = "homer2", homerfile = homerbin, genome = "error"))
        expect_is(findMotifHits(pwmL, unname(gr), min.score = "90%", method = "homer2", homerfile = homerbin, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10), "GRanges")

        expect_error(findMotifHits(pwmL, sf,   min.score = TRUE, method = "homer2", homerfile = homerbin))
        expect_error(findMotifHits(pwmL, seqs, min.score = TRUE, method = "homer2", homerfile = homerbin))

        res1h  <- findMotifHits(tf,   sf,   min.score = "90%", method = "homer2", homerfile = homerbin) # character,character
        res2h  <- findMotifHits(tf,   seq1, min.score = "90%", method = "homer2", homerfile = homerbin) # character,DNAString
        res3h  <- findMotifHits(tf,   seqs, min.score = "90%", method = "homer2", homerfile = homerbin) # character,DNAStringSet
        res4h  <- findMotifHits(pwm,  sf,   min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrix,character
        res5h  <- findMotifHits(pwm,  seq1, min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrix,DNAString
        res6h  <- findMotifHits(pwm,  seqs, min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrix,DNAStringSet
        res7h  <- findMotifHits(pwmL, sf,   min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrixList,character
        res7hb <- findMotifHits(pwmL, sf,   min.score = 4.456, method = "homer2", homerfile = homerbin) # PWMatrixList,character
        res8h  <- findMotifHits(pwmL, seq1, min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrixList,DNAString
        res9h  <- findMotifHits(pwmL, seqs, min.score = "90%", method = "homer2", homerfile = homerbin) # PWMatrixList,DNAStringSet
        res9hb <- findMotifHits(pwmL, seqs, min.score = 4.456, method = "homer2", homerfile = homerbin) # PWMatrixList,DNAStringSet
        res10h <- findMotifHits(pwm,  gr,   min.score = "90%", method = "homer2", homerfile = homerbin, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrix,GRanges
        res11h <- findMotifHits(pwmL, gr,   min.score = "90%", method = "homer2", homerfile = homerbin, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) # PWMatrixList,GRanges

        # consistency between "matchPWM" and "homer2"
        expect_equivalent(res1, res1h)
        expect_equivalent(res2, res2h)

        # consistency between methods
        expect_equal(res1h, res3h)
        expect_equal(res1h, res4h)
        expect_equal(res2h, res5h)
        expect_equal(res1h, res6h)
        expect_equal(res1h, res7h)
        expect_equal(res1h, res7hb)
        expect_equal(res2h, res8h)
        expect_equal(res1h, res9h)
        expect_equal(res1h, res9hb)
        expect_equal(res1h, res10h)
        expect_equal(res1h, res11h)
    }

    unlink(c(sf, tf))
})

