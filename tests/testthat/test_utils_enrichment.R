context("utils_enrichment")

test_that(".binomEnrichmentTests works", {
    res <- .binomEnrichmentTest(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, verbose = TRUE)
    expect_equal(res[1], 
                 log(binom.test(x = 9, n = 20, p = 5/17, 
                                alternative = "greater")$p.value))
    expect_equal(res[2], 
                 log(binom.test(x = 3, n = 20, p = 8/17, 
                                alternative = "greater")$p.value))
})

test_that(".fisherEnrichmentTest works", {
    res <- .fisherEnrichmentTest(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, verbose = TRUE)
    expect_equal(res[1], 
                 log(fisher.test(x = matrix(c(9, 20 - 9, 5, 17 - 5), nrow = 2),  
                                 alternative = "greater")$p.value))
    expect_equal(res[2], 
                 log(fisher.test(x = matrix(c(3, 20 - 3, 8, 17 - 8), nrow = 2), 
                                 alternative = "greater")$p.value))
})

test_that(".calcPearsonResiduals works", {
    res <- .calcPearsonResiduals(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20)
    expect_equal(res[1], 
                 (9 - (9 + 5)/(20 + 17) * 20)/sqrt((9 + 5)/(20 + 17) * 20 * 
                                                       (1 - 20/(20 + 17)) * 
                                                       (1 - (5 + 9)/(20 + 17))))
    expect_equal(res[2], 
                 (3 - (3 + 8)/(20 + 17) * 20)/sqrt((3 + 8)/(20 + 17) * 20 * 
                                                       (1 - 20/(20 + 17)) * 
                                                       (1 - (8 + 3)/(20 + 17))))
    
})

test_that(".calcExpFg works", {
    res <- .calcExpFg(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20)
    expect_equal(res[1], (9 + 5)/(20 + 17) * 20)
    expect_equal(res[2], (3 + 8)/(20 + 17) * 20)
})

test_that(".calcLog2Enr works", {
    res <- .calcLog2Enr(
        matchCountBg = c(5, 8), totalWeightBg = 17, 
        matchCountFg = c(9, 3), totalWeightFg = 20, pseudocount = 3)
    expect_equal(res[1], log2((9/20 * 17 + 3)/(5/17 * 17 + 3)))
    expect_equal(res[2], log2((3/20 * 17 + 3)/(8/17 * 17 + 3)))
})


