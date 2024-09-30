context("plotting")

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

test_that("plotBinDiagnostics() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)
    
    expect_is(plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "length"), 
              "list")
    expect_is(plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "GCfrac"), 
              "list")
    expect_is(plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "dinucfreq"), 
              "Heatmap")
    
    dev.off()
    unlink(tf)
    
    expect_error(plotBinDiagnostics(seqs = x, bins = b1))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = as.numeric(b1)))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = as.character(b1)))
    expect_error(plotBinDiagnostics(seqs = seqs, bins = b1, aspect = "missing"))
})

test_that("plotBinScatter() runs", {
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_true(plotBinScatter(x = x, y = x, b = b1))
    expect_error(plotBinScatter(x = x, y = x, b = b1, cols = "gray"))
    expect_true(plotBinScatter(x = x, y = x, b = b1, cols = "gray", legend = FALSE))

    dev.off()
    unlink(tf)
})


test_that("plotMotifHeatmaps() runs", {
    expect_error(plotMotifHeatmaps(x = se, cluster = "error"))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_is(plotMotifHeatmaps(x = se, which.plots = "pearsonResid", cluster = FALSE, show_motif_GC = TRUE), "list")
    expect_is(plotMotifHeatmaps(x = se, which.plots = "negLog10Padj", cluster = TRUE, show_seqlogo = TRUE), "list")
    cl <- hclust(dist(SummarizedExperiment::assay(se, "log2enr")))
    expect_is(plotMotifHeatmaps(x = se, which.plots = "log2enr", cluster = cl, show_dendrogram = TRUE), "list")

    se2 <- se
    tmp <- SummarizedExperiment::assay(se2, "pearsonResid")
    tmp[1:2, ] <- NA
    SummarizedExperiment::assay(se2, "pearsonResid") <- tmp
    expect_warning(res <- plotMotifHeatmaps(x = se2, which.plots = "log2enr", cluster = TRUE))
    expect_is(res, "list")
    
    expect_error(plotMotifHeatmaps(x = se, show_bin_legend = "error"))
    
    result_true <- plotMotifHeatmaps(x = se, which.plots = "pearsonResid", show_bin_legend = TRUE, doPlot = FALSE)
    expect_true(result_true$pearsonResid@top_annotation@anno_list$bin@show_legend)
    
    result_false <- plotMotifHeatmaps(x = se, which.plots = "pearsonResid", show_bin_legend = FALSE, doPlot = FALSE)
    expect_false(result_false$pearsonResid@top_annotation@anno_list$bin@show_legend)
    
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

    pdf(file = tf)

    expect_error(plotSelectionProb(se = "error", selProbMin = 0.5))
    expect_error(plotSelectionProb(se = ss, directional = "error"), "logical")
    expect_error(plotSelectionProb(se = ss, directional = TRUE, selProbMin = 2.0), "within")
    expect_error(plotSelectionProb(se = ss, selProbMinPlot = "error"), "numeric")
    expect_error(plotSelectionProb(se = ss, selProbMin = 0.5, selProbMinPlot = 0.6))
    expect_error(plotSelectionProb(se = ss, showSelProbMin = "error"))
    expect_error(plotSelectionProb(se = ss, col = "error"), "length 3")
    expect_error(plotSelectionProb(se = ss, method = "error"), "should be one of")
    expect_error(plotSelectionProb(se = ss, legend = "error"))

    expect_null(plotSelectionProb(ss, selProbMin = 1.0, selProbMinPlot = 0.99))
    expect_is(plotSelectionProb(ss), "matrix")
    expect_is(plotSelectionProb(ss, FALSE), "matrix")
    
    dev.off()
    unlink(tf)
})

