context("binning")

test_that("bin() works properly", {
    # normal distribution
    set.seed(1)
    x <- rnorm(1000)

    b1 <- bin(x, binmode = "equalN", nElements = 100)
    expect_equal(nlevels(b1), 10L)
    expect_true(all(table(b1) == 100))
    expect_equal(getZeroBin(b1), NA)

    b2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
    expect_true(table(b2)[attr(b2, "bin0")] == 400L)
    expect_equal(getZeroBin(b2), 7)
    b2p <- setZeroBin(b2, 3)
    expect_true(table(b2p)[attr(b2p, "bin0")] == 50L)
    expect_equal(getZeroBin(b2p), 3)

    b3 <- bin(x, binmode = "equalWidth", nBins = 5)
    expect_equal(nlevels(b3), 5L)
    expect_equal(diff(attr(b3, "breaks")), rep(1.363665, 5L), tolerance = 1e-6)

    b4 <- bin(x, binmode = "equalWidth", nBins = 5, minAbsX = 0.6)
    expect_true(table(b4)[attr(b4, "bin0")] == 440L)
    
    # asymmetric distribution
    set.seed(2)
    x <- c(rnorm(225, 0.006424305, 0.0792525),
           rnorm(775, 0.735741802, 0.5800667))
    
    a1 <- bin(x, binmode = "equalN", nElements = 100)
    expect_equal(nlevels(a1), 10L)
    expect_true(all(table(a1) == 100))
    
    a2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
    expect_true(table(a2)[attr(a2, "bin0")] == 550L)
    
    a3 <- bin(x, binmode = "equalWidth", nBins = 5)
    expect_equal(nlevels(a3), 5L)
    expect_equal(diff(attr(a3, "breaks")), rep(0.6648299, 5L), tolerance = 1e-6)
    
    a4 <- bin(x, binmode = "equalWidth", nBins = 5, minAbsX = 0.6)
    expect_true(table(a4)[attr(a4, "bin0")] == 526L)
})

