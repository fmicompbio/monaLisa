context("region sampling")

test_that("sample_random_regions() works properly", {
    gr <- GenomicRanges::GRanges(seqnames = paste0("chr", seq.int(20L)),
                                 ranges = IRanges::IRanges(start = 1L, width = 10000L))

    expect_error(sample_random_regions())
    expect_error(sample_random_regions(N = 1:2))
    expect_error(sample_random_regions(N = "error"))
    # todo: add some more

    res1 <- sample_random_regions(mappableRegions = gr)
    res2 <- sample_random_regions(mappableRegions = gr)

    expect_is(res1, "GRanges")
    expect_identical(res1, res2)
    expect_length(res1, 100L)
    expect_true(all(width(res1) == 200L))
    # todo: add some more (e.g. for fractionCGI)
})

