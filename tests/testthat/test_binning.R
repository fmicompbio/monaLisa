test_that("bin() works properly", {
    set.seed(1)
    x <- rnorm(1000)
    b1 <- bin(x, binmode = "equalN", nElements = 100)
    expect_equal(nlevels(b1), 10L)
    expect_true(all(table(b1) == 100))
    b2 <- bin(x, binmode = "equalWidth", nBins = 5)
    expect_equal(nlevels(b2), 5L)
    expect_equal(diff(attr(b2, "breaks")), rep(1.363665, 5L), tolerance = 1e-6)
    b3 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
    expect_true(table(b3)[attr(b3, "bin0")] == 400L)
})

