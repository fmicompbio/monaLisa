context("motifs")

test_that("motifs can be written to/read from files", {
    # create dummy motif file and PWMatrixList
    m <- matrix(data = c(.1,.3,.3,.3,.01,.01,.01,.97,.5,.0,.0,.5), nrow = 4, dimnames = list(c("A","C","G","T"), NULL))
    pwm <- TFBSTools::PWMatrixList(TFBSTools::PWMatrix(ID = "CTA", name = "name", profileMatrix = log2((m+1e-9)/(.25+1e-9))))
    pwm <- pwm[c(1,1,1)]
    names(pwm) <- rep("name",3)
    tmpf1 <- tempfile()
    tmpf2 <- tempfile()
    tmpstr <- ">CTA\tname\t4.80\n0.1\t0.3\t0.3\t0.3\n0.01\t0.01\t0.01\t0.97\n0.5\t0\t0\t0.5\n"
    cat(tmpstr[c(1,1,1)], file = tmpf1, sep = "")

    # read from file
    pwm.file <- monaLisa:::.readPWMsFromHomer2File(tmpf1)
    expect_equal(pwm, pwm.file)

    # write to file
    monaLisa:::.dumpPWMsToHomer2File(pwmL = pwm, fname = tmpf2, absscore = 6.93)
    lns <- readLines(tmpf2)
    expect_identical(">CTA\tname\t4.803510", lns[1])
    m.file <- as.matrix(read.delim(tmpf2, header = FALSE, nrows = 3, skip = 1))
    expect_equivalent(m, t(m.file))

    # clean up
    unlink(x = c(tmpf1, tmpf2))
})

