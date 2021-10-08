context("region sampling")

test_that("sampleRandomRegions() works properly", {
    gr <- GenomicRanges::GRanges(seqnames = paste0("chr", seq.int(20L)),
                                 ranges = IRanges::IRanges(start = 1L, width = 10000L))

    expect_error(sampleRandomRegions())
    expect_error(sampleRandomRegions(N = 1:2))
    expect_error(sampleRandomRegions(N = "error"))
    # todo: add some more

    set.seed(123)
    res1 <- sampleRandomRegions(allowedRegions = gr)
    set.seed(123)
    res2 <- sampleRandomRegions(allowedRegions = gr)

    expect_is(res1, "GRanges")
    expect_identical(res1, res2)
    expect_length(res1, 100L)
    expect_true(all(width(res1) == 200L))
    # todo: add some more (e.g. for fractionCGI)
})

