test_that("randLassoStabSel() works properly", {

    # create data set
    set.seed(555)
    Y <- rnorm(n = 500, mean = 2, sd = 1)
    X <- matrix(data = NA, nrow = length(Y), ncol = 50)
    for (i in seq_len(ncol(X))) {
        X[ ,i] <- runif(n = 500, min = 0, max = 3)
    }
    s_cols <- sample(x = seq_len(ncol(X)), size = 10, replace = FALSE)
    for (i in seq_along(s_cols)) {
        X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
    }

    X2 <- X
    Y2 <- Y
    rownames(X2) <- paste0("peak", seq_len(nrow(X2)))
    colnames(X2) <- paste0("motif", seq_len(ncol(X2)))
    names(Y2) <- paste0("peak", seq_along(Y2))

    # randomized lasso stability selection
    ss <- monaLisa::randLassoStabSel(x = X, y = Y)

    # tests
    expect_s4_class(ss, "SummarizedExperiment")
    expect_identical(rowData(ss)$y, Y)
    expect_identical(ss$selProb, colData(ss)[, ncol(colData(ss))])
    expect_identical(ss$selAUC,
                     rowMeans(as.matrix(colData(ss)[, grep("^regStep", colnames(colData(ss)))])))
    expect_true(all(ss$selProb >= 0 & ss$selProb <= 1))
    expect_true(all(ss$selAUC >= 0 & ss$selAUC <= 1))
    expect_true(all(s_cols %in% metadata(ss)$stabsel.params.selected))
    expect_identical(dim(ss), c(500L, 50L))
    expect_identical(length(Y), nrow(ss))
    expect_true(!is.null(SummarizedExperiment::assay(ss)))
    expect_error(randLassoStabSel(x = as.data.frame(X), y = Y))
    expect_error(randLassoStabSel(x = X2[1:100, ], y = Y2[2:101]))
    expect_error(randLassoStabSel(x = X2[1:100, ], y = Y2[2:100]))
    
    # with a different family
    set.seed(123)
    ssbnf <- randLassoStabSel(x = X, y = as.numeric(Y > 2), 
                              glmnet.args = list(family = "binomial"))
    expect_s4_class(ssbnf, "SummarizedExperiment")
    expect_identical(rowData(ssbnf)$y, as.numeric(Y > 2))
    expect_identical(ssbnf$selProb, colData(ssbnf)[, ncol(colData(ssbnf))])
    expect_identical(ssbnf$selAUC,
                     rowMeans(as.matrix(colData(ssbnf)[, grep("^regStep", colnames(colData(ssbnf)))])))
    expect_true(all(ssbnf$selProb >= 0 & ssbnf$selProb <= 1))
    expect_true(all(ssbnf$selAUC >= 0 & ssbnf$selAUC <= 1))
    expect_true(all(metadata(ssbnf)$stabsel.params.selected %in% s_cols))
    expect_identical(dim(ssbnf), c(500L, 50L))
    expect_identical(length(Y), nrow(ssbnf))
    expect_true(!is.null(SummarizedExperiment::assay(ssbnf)))
    
    # ignore selected arguments if passed to glmnet.args
    set.seed(123)
    expect_warning({
        ssbnf2 <- randLassoStabSel(x = X, y = as.numeric(Y > 2), 
                                   glmnet.args = list(family = "binomial", 
                                                      x = 3, weakness = 0.1))
    }, "Ignoring the following elements")
    expect_identical(ssbnf, ssbnf2)
})

test_that("randLassoStabSel() is deterministic", {

  # create data set
  set.seed(555)
  Y <- rnorm(n = 500, mean = 2, sd = 1)
  X <- matrix(data = NA, nrow = length(Y), ncol = 50)
  for (i in seq_len(ncol(X))) {
    X[ ,i] <- runif(n = 500, min = 0, max = 3)
  }
  s_cols <- sample(x = seq_len(ncol(X)), size = 10, replace = FALSE)
  for (i in seq_along(s_cols)) {
    X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
  }

  # randomized lasso stability selection
  set.seed(123)
  ss1 <- monaLisa::randLassoStabSel(x = X, y = Y)
  set.seed(123)
  ss2 <- monaLisa::randLassoStabSel(x = X, y = Y)

  # tests
  expect_identical(ss1, ss2)

})


test_that(".glmnetRandomizedLasso() works properly", {

    # create data set
    set.seed(555)
    Y <- rnorm(n = 500, mean = 2, sd = 1)
    X <- matrix(data = NA, nrow = length(Y), ncol = 50)
    for (i in seq_len(ncol(X))) {
        X[ ,i] <- runif(n = 500, min = 0, max = 3)
    }
    s_cols <- sample(x = ncol(X), size = 10, replace = FALSE)
    for (i in seq_along(s_cols)) {
        X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
    }

    # tests
    # ... x as data.frame
    expect_warning(
      expect_message(.glmnetRandomizedLasso(x = as.data.frame(X), y = Y, q = 11),
                     "coerced to a model matrix without intercept"),
      "Number of nonzero coefficients along the path exceeds")

    # ... with specific lambda
    expect_error(.glmnetRandomizedLasso(x = X, y = Y, q = 11, lambda = 5))

    # ... with type="anticonservative"
    set.seed(123)
    rl <- .glmnetRandomizedLasso(x = X, y = Y, q = 11, type = "anticonservative")
    expect_type(rl, "list")
    expect_identical(names(rl), c("selected", "path"))
    expect_type(rl$selected, "logical")
    expect_type(rl$path, "logical")

    # expected number of selected variables
    expect_identical(sum(rl$selected), 11L)
    
    # with different family
    set.seed(123)
    expect_warning(
        expect_message(
            bnf <- .glmnetRandomizedLasso(x = as.data.frame(X), 
                                          y = as.numeric(Y > 0), 
                                          q = 11, family = "binomial"),
            "coerced to a model matrix without intercept"),
        "Number of nonzero coefficients along the path exceeds")
    expect_identical(sum(bnf$selected), 10L)

})

