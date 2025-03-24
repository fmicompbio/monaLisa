test_that("bin() works properly", {
    # normal distribution
    set.seed(1)
    x <- rnorm(1000)

    b1 <- bin(x, binmode = "equalN", nElements = 100)
    expect_equal(nlevels(b1), 10L)
    expect_true(all(table(b1) == 100))
    expect_equal(getZeroBin(b1), NA)

    b2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
    expect_true(table(b2)[getZeroBin(b2)] == 400L)
    expect_equal(getZeroBin(b2), 7)
    b2p <- setZeroBin(b2, 3)
    expect_true(table(b2p)[getZeroBin(b2p)] == 50L)
    expect_equal(getZeroBin(b2p), 3)

    b3 <- bin(x, binmode = "equalWidth", nBins = 5)
    expect_equal(nlevels(b3), 5L)
    expect_equal(diff(attr(b3, "breaks")), rep(1.363665, 5L), tolerance = 1e-6)

    b4 <- bin(x, binmode = "equalWidth", nBins = 5, minAbsX = 0.6)
    expect_true(table(b4)[getZeroBin(b4)] == 440L)

    expect_warning(b5 <- bin(x, binmode = "equalN", nBins = 5, minAbsX = 3))

    b6 <- bin(x = seq(-20, 10), binmode = "equalWidth", nBins = 5, minAbsX = 10)
    expect_equal(nlevels(b6), 5L)
    expect_equal(diff(attr(b6, "breaks")), c(rep(2.5, 4L), 20), tolerance = 1e-6)
    
    # uniform distribution, all negative values
    u1 <- bin(seq(-20, -10), binmode = "equalWidth", nBins = 5)
    expect_type(u1, "integer")
    expect_true(all(c("binmode", "breaks", "bin0") %in% names(attributes(u1))))
    
    # uniform distribution, all negative values, set minAbsX
    expect_warning({
        u2 <- bin(seq(-20, -10), binmode = "equalWidth", nBins = 5, minAbsX = 5)
    }, "All data values are outside")
    expect_type(u2, "integer")
    expect_true(all(c("binmode", "breaks", "bin0") %in% names(attributes(u2))))
    u3 <- bin(seq(-20, -10), binmode = "equalWidth", nBins = 5, minAbsX = NULL)
    expect_identical(u2, u3)

    # uniform distribution, all negative values, set minAbsX (equalN)
    expect_warning({
        u2 <- bin(seq(-20, -10), binmode = "equalN", nBins = 5, minAbsX = 5)
    }, "All data values are outside")
    expect_type(u2, "integer")
    expect_true(all(c("binmode", "breaks", "bin0") %in% names(attributes(u2))))
    u3 <- bin(seq(-20, -10), binmode = "equalN", nBins = 5, minAbsX = NULL)
    expect_identical(u2, u3)
    
    # asymmetric distribution
    set.seed(2)
    x <- c(rnorm(225, 0.006424305, 0.0792525),
           rnorm(775, 0.735741802, 0.5800667))

    a1 <- bin(x, binmode = "equalN", nElements = 100)
    expect_equal(nlevels(a1), 10L)
    expect_true(all(table(a1) == 100))

    expect_warning(a2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6))
    expect_true(table(a2)[getZeroBin(a2)] == 550L)

    a3 <- bin(x, binmode = "equalWidth", nBins = 5)
    expect_equal(nlevels(a3), 5L)
    expect_equal(diff(attr(a3, "breaks")), rep(0.6648299, 5L), tolerance = 1e-6)

    expect_warning(a4 <- bin(x, binmode = "equalWidth", nBins = 5, minAbsX = 0.6))
    expect_true(table(a4)[getZeroBin(a4)] == 526L)
})

test_that("getZeroBin() and setZeroBin() work properly", {
    set.seed(1)
    x <- rnorm(1000)

    b1 <- bin(x, nElements = 100)
    b2 <- bin(x, nElements = 100, minAbsX = 0.5)
    b3 <- cut(x, breaks = 10)

    expect_true(is.na(getZeroBin(b1)))
    expect_identical(getZeroBin(b2), 4L)
    expect_true(is.null(getZeroBin(b3)))

    expect_error(setZeroBin(b1, FALSE))
    expect_error(setZeroBin(b1, 2000))
    expect_error(setZeroBin(b1, "error"))

    expect_type(b4 <- setZeroBin(b1, 6), "integer")
    expect_type(b5 <- setZeroBin(b1, "(-0.0353,0.245]"), "integer")
    expect_type(b6 <- setZeroBin(b1, NA), "integer")

    expect_type(getZeroBin(b4), "integer")
    expect_identical(getZeroBin(b4), getZeroBin(b5))
    expect_true(is.na(getZeroBin(b6)))
})

