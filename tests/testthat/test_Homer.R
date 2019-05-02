context("Homer")

test_that("findHomer() works properly", {
    res <- findHomer("I-do-not-exist")
    expect_true(is.na(res))

    res <- findHomer("resL.rds", dirs = system.file("extdata", package = "lisa"))
    expect_true(file.exists(res))
})

test_that("dumpJaspar() works properly", {
    tmp1 <- tempfile()

    expect_error(dumpJaspar(filename = system.file("extdata", "resL.rds", package = "lisa")))
    expect_error(dumpJaspar(filename = tmp1, relScoreCutoff = "error"))
    expect_true(dumpJaspar(filename = tmp1, opts = list(ID = c("MA0006.1", "MA0007.3"))))

    unlink(tmp1)
})

test_that("prepareHomer() works properly", {
    gr <- GenomicRanges::GRanges(rep(c("chr1","chr2","chr3"), each = 10),
                                 IRanges::IRanges(start = rep(10 * 1:10, 3), width = 5L))
    b <- 1:30 %% 3 + 1
    bF <- factor(b)
    tmp1 <- tempfile()
    dir.create(tmp1)
    tmp2 <- tempfile()
    fname <- system.file("extdata", "resL.rds", package = "lisa")

    expect_error(prepareHomer(gr = gr, b = b, genomedir = "genomedir", outdir = tmp1,
                              motifFile = fname, homerfile = fname, regionsize = "given", Ncpu = 2))
    expect_error(prepareHomer(gr = as(gr, "data.frame"), b = bF[1:10], genomedir = "genomedir", outdir = tmp2,
                              motifFile = fname, homerfile = fname, regionsize = "given", Ncpu = 2))
    expect_error(prepareHomer(gr = gr, b = bF, genomedir = "genomedir", outdir = tmp2,
                              motifFile = "error", homerfile = fname, regionsize = "given", Ncpu = 2))

    expect_identical(prepareHomer(gr = gr, b = bF, genomedir = "genomedir", outdir = tmp2,
                                  motifFile = fname, homerfile = fname, regionsize = "given", Ncpu = 2),
                     file.path(tmp2, "run.sh"))

    unlink(c(tmp1, tmp2), recursive = TRUE, force = TRUE)
})

test_that("parseHomerOutput() works properly", {
    outfile <- system.file("extdata", "homer_output.txt.gz", package = "lisa")

    res <- parseHomerOutput(c(outfile, outfile))
    expect_length(res, 4L)
    expect_identical(names(res), c("p", "FDR", "enr", "log2enr"))
    expect_identical(res$p[,1], res$p[,2])
    expect_true(all(sapply(res, dim) == c(519L, 2L)))
    expect_equal(sum(res$enr), -914.6696)
})
