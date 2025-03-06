# create data

# ... binning
set.seed(1)
x <- rnorm(1000)
b1 <- bin(x, binmode = "equalN", nElements = 100)
b2 <- bin(x, binmode = "equalN", nElements = 50, minAbsX = 0.6)
se <- readRDS(system.file("extdata", "results.binned_motif_enrichment_LMRs.rds",
                          package = "monaLisa"))[1:10, 1:8]
seqs <- Biostrings::DNAStringSet(
    vapply(seq_along(x),
           function(i) paste(sample(c("A", "C", "G", "T"), 10,
                                    replace = TRUE), collapse = ""), "")
)

# ... stability selection
Y <- rnorm(n = 100, mean = 2, sd = 1)
X <- matrix(data = runif(n = 20 * 100, min = 0, max = 3), nrow = length(Y), ncol = 20)
for (i in sample(x = seq_len(ncol(X)), size = 10, replace = FALSE))
    X[ ,i] <- X[ ,i] + Y * c(1, -1)[(i %% 2) + 1]
ss <- monaLisa::randLassoStabSel(x = X, y = Y)


test_that("getColsByBin() works properly", {
    c1 <- getColsByBin(as.numeric(b1))
    c2 <- getColsByBin(b2)

    expect_length(c1, 1000L)
    expect_identical(as.vector(c1), as.vector(getColsByBin(c1)))
    expect_equal(sort(unname(table(c2))), sort(unname(table(b2))))
})


test_that("plotBinHist() runs", {
    expect_warning(plotBinHist(x = x, b = b1, legend = "topright"))
    expect_warning(plotBinHist(x = x, b = b1, legend.cex = 1.0))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_s3_class(plotBinHist(x = x, b = b1), "ggplot")

    dev.off()
    unlink(tf)
})


test_that("plotBinDensity() runs", {
    expect_warning(plotBinDensity(x = x, b = b1, legend = "topright"))
    expect_warning(plotBinDensity(x = x, b = b1, legend.cex = 1.0))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_s3_class(plotBinDensity(x = x, b = b1), "ggplot")

    dev.off()
    unlink(tf)
})

test_that("plotBinDiagnostics() runs", {
    expect_error(plotBinDiagnostics(seqs = seqs, bins = b1[1:3]))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_s3_class(
        plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "length"),
        "ggplot")
    expect_s3_class(
        plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "GCfrac"),
        "ggplot")
    expect_s4_class(
        plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "dinucfreq"),
        "Heatmap")

    dev.off()
    unlink(tf)

    expect_error(plotBinDiagnostics(seqs = x, bins = b1))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = as.numeric(b1)))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = as.character(b1)))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "missing"))
})

test_that("plotBinScatter() runs", {
    expect_warning(plotBinScatter(x = x, y = x, b = b1, cols = "red",
                                  legendPosition = "right"))
    expect_warning(plotBinScatter(x = x, y = x, b = b1, legend = "topright"))
    expect_warning(plotBinScatter(x = x, y = x, b = b1, legend.cex = 1.0))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_s3_class(plotBinScatter(x = x, y = x, b = b1), "ggplot")
    expect_s3_class(plotBinScatter(x = x, y = x, b = b1,
                                   cols = "gray",
                                   legendPosition = "none"), "ggplot")

    dev.off()
    unlink(tf)
})


test_that("plotMotifHeatmaps() runs", {
    expect_error(plotMotifHeatmaps(x = se, cluster = "error"))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_type(plotMotifHeatmaps(x = se, which.plots = "pearsonResid", cluster = FALSE, show_motif_GC = TRUE), "list")
    expect_type(plotMotifHeatmaps(x = se, which.plots = "negLog10Padj", cluster = TRUE, show_seqlogo = TRUE), "list")
    cl <- hclust(dist(SummarizedExperiment::assay(se, "log2enr")))
    expect_type(plotMotifHeatmaps(x = se, which.plots = "log2enr", cluster = cl, show_dendrogram = TRUE), "list")

    se2 <- se
    tmp <- SummarizedExperiment::assay(se2, "pearsonResid")
    tmp[1:2, ] <- NA
    SummarizedExperiment::assay(se2, "pearsonResid") <- tmp
    expect_warning(res <- plotMotifHeatmaps(x = se2, which.plots = "log2enr", cluster = TRUE))
    expect_type(res, "list")

    expect_error(plotMotifHeatmaps(x = se, show_bin_legend = "error"))

    result_true <- plotMotifHeatmaps(x = se, which.plots = "pearsonResid",
                                     show_bin_legend = TRUE, doPlot = FALSE,
                                     highlight = rep(c(TRUE, FALSE), c(3, 7)))
    expect_true(result_true$pearsonResid@top_annotation@anno_list$bin@show_legend)

    result_false <- plotMotifHeatmaps(x = se, which.plots = c("pearsonResid", "negLog10P"),
                                      show_bin_legend = FALSE, doPlot = FALSE,
                                      maxEnr = 4, maxSig = 8)
    expect_false(result_false$pearsonResid@top_annotation@anno_list$bin@show_legend)

    dev.off()
    unlink(tf)
})


test_that("plotStabilityPaths() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_error(plotStabilityPaths("error"))
    expect_s3_class(plotStabilityPaths(ss), "ggplot")

    dev.off()
    unlink(tf)
})


test_that("plotSelectionProb() runs", {
    tf <- tempfile(fileext = ".pdf")

    pdf(file = tf)

    expect_error(plotSelectionProb(se = "error", selProbMin = 0.5))
    expect_error(plotSelectionProb(se = ss, directional = "error"), "logical")
    expect_error(plotSelectionProb(se = ss, directional = TRUE, selProbMin = 2.0), "between")
    expect_error(plotSelectionProb(se = ss, selProbMinPlot = "error"), "numeric")
    expect_error(plotSelectionProb(se = ss, selProbMin = 0.5, selProbMinPlot = 0.6))
    expect_error(plotSelectionProb(se = ss, showSelProbMin = "error"), "logical")
    expect_error(plotSelectionProb(se = ss, selColor = "error"))
    expect_error(plotSelectionProb(se = ss, method = "error"), "should be one of")
    expect_s3_class(plotSelectionProb(ss), "ggplot")

    dev.off()
    unlink(tf)
})

