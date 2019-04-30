context("Homer")

test_that("findHomer() works properly", {
    res <- findHomer("NAMESPACE", dirs = system.file(package = "lisa"))
    expect_true(file.exists(res))
})

