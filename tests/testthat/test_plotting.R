context("plotting")

# create data

# ... binning
set.seed(1)
x <- rnorm(1000)
b1 <- bin(x, binmode = "equalN", nElements = 100)
b2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
se <- readRDS(system.file("extdata", "se.rds", package = "monaLisa"))[1:10, 1:8]

# ... stability selection
Y <- rnorm(n = 100, mean = 2, sd = 1)
X <- matrix(data = runif(n = 20 * 100, min = 0, max = 3), nrow = length(Y), ncol = 20)
for (i in sample(x = 1:ncol(X), size = 10, replace = FALSE))
    X[ ,i] <- X[ ,i] + Y
ss <- monaLisa::randomized_stabsel(x = X, y = Y)


test_that("getColsByBin() works properly", {
    c1 <- getColsByBin(as.numeric(b1))
    c2 <- getColsByBin(b2)

    expect_length(c1, 1000L)
    expect_identical(as.vector(c1), as.vector(getColsByBin(c1)))
    expect_equal(sort(unname(table(c2))), sort(unname(table(b2))))
})


test_that("plotBinHist() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_is(plotBinHist(x = x, b = b1), "histogram")

    dev.off()
    unlink(tf)
})


test_that("plotBinDensity() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_is(plotBinDensity(x = x, b = b1), "density")

    dev.off()
    unlink(tf)
})


test_that("plotBinScatter() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(plotBinScatter(x = x, y = x, b = b1))
    expect_error(plotBinScatter(x = x, y = x, b = b1, cols = "gray"))
    expect_null(plotBinScatter(x = x, y = x, b = b1, cols = "gray", legend = FALSE))

    dev.off()
    unlink(tf)
})


test_that("plotMotifHeatmaps() runs", {
    expect_error(plotMotifHeatmaps(x = se, cluster = "error"))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_is(plotMotifHeatmaps(x = se, which.plots = "enr", cluster = FALSE, show_motif_GC = TRUE), "list")
    expect_is(plotMotifHeatmaps(x = se, which.plots = "FDR", cluster = TRUE, show_seqlogo = TRUE), "list")
    cl <- hclust(dist(SummarizedExperiment::assay(se, "log2enr")))
    expect_is(plotMotifHeatmaps(x = se, which.plots = "log2enr", cluster = cl, show_dendrogram = TRUE), "list")

    dev.off()
    unlink(tf)
})


test_that("plotStabilityPaths() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_error(plotStabilityPaths("error"))
    expect_true(plotStabilityPaths(ss))

    dev.off()
    unlink(tf)
})


test_that("plotSelectionProb() runs", {
    tf <- tempfile(fileext = ".pdf")

    ss2 <- ss
    # SummarizedExperiment::rowData(ss2)[,1] <- FALSE

    pdf(file = tf)

    expect_error(plotSelectionProb("error"))
    # expect_error(plotSelectionProb(ss2, onlySelected = TRUE))
    expect_true(plotSelectionProb(ss, onlySelected = FALSE))
    expect_true(plotSelectionProb(ss, onlySelected = TRUE))

    dev.off()
    unlink(tf)
})


test_that("plotMotifDirectionality() runs", {
    tf <- tempfile(fileext = ".pdf")
    
    X2 <- X
    colnames(X2) <- paste0("pred", 1:ncol(X2))
    
    pdf(file = tf)
    
    expect_error(plotMotifDirectionality("error"))
    expect_error(plotMotifDirectionality(se = NULL))
    expect_true(plotMotifDirectionality(se = ss))
    
    dev.off()
    unlink(tf)
})


