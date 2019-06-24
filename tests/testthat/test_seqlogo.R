context("seqlogo")

m <- matrix(data = c(.1,.3,.3,.3,.01,.01,.01,.97,.5,.0,.0,.5), nrow = 4,
            dimnames = list(c("A","C","G","T"), NULL))
pfm <- TFBSTools::PFMatrix(ID = "test", name = "name", profileMatrix = 100 * m)

test_that("letter adding functions work as expected", {
    rA <- letterA(x.pos = 0, y.pos = 0, ht = 1, wt = 1)
    rC <- letterC(x.pos = 0, y.pos = 0, ht = 1, wt = 1)
    rG <- letterG(x.pos = 0, y.pos = 0, ht = 1, wt = 1)
    rT <- letterT(x.pos = 0, y.pos = 0, ht = 1, wt = 1)

    expect_is(rA, "list")
    expect_is(rC, "list")
    expect_is(rG, "list")
    expect_is(rT, "list")

    expect_identical(names(rA), c("x", "y", "id", "fill"))
    expect_identical(names(rA), names(rC))
    expect_identical(names(rA), names(rG))
    expect_identical(names(rA), names(rT))

    expect_identical(unname(lengths(rA)), c( 13L,  13L,  13L, 2L))
    expect_identical(unname(lengths(rC)), c(778L, 778L, 778L, 1L))
    expect_identical(unname(lengths(rG)), c(785L, 785L, 785L, 2L))
    expect_identical(unname(lengths(rT)), c(  8L,   8L,   8L, 1L))

    expect_identical(addLetter(letters = list(), which = "A", x.pos = 0, y.pos = 0, ht = 1, wt = 1), rA)
})

test_that("information content in a  weight matrix is calculated correctly", {
    expect_identical(pfm2ic(m), TFBSTools::totalIC(TFBSTools::toICM(pfm, pseudocounts = 0)))
})

test_that("sequence logo can be drawn", {
    g1 <- seqLogoGrob(pfm, xmax = 5L)
    g2 <- seqLogoGrob(pfm, xmax = 5L, xjust = "center")
    g3 <- seqLogoGrob(pfm, xmax = 5L, xjust = "right")
    hmr <- anno_seqlogo(list(g1, g2, g3), which = "row")
    hmc <- anno_seqlogo(list(g1, g2, g3), which = "column")

    expect_is(g1, "polygon")
    expect_is(hmr, "AnnotationFunction")
    expect_is(hmc, "AnnotationFunction")

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(ComplexHeatmap::draw(hmr, test = "Test"))
    expect_null(ComplexHeatmap::draw(hmc, test = "Test"))

    dev.off()
    unlink(tf)
})

