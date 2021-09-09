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
    expect_true(is(ss, "SummarizedExperiment"))
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
    rl <- .glmnetRandomizedLasso(x = X, y = Y, q = 11, type = "anticonservative")
    expect_is(rl, "list")
    expect_identical(names(rl), c("selected", "path"))
    expect_is(rl$selected, "logical")
    expect_is(rl$path, "matrix")
  
    # outputs differ depending on R version due to random number generator with 'sample' function (check RNGkind) --> different as of R 3.6.0
    Rmajor <- as.numeric(R.version$major)
    Rminor <- as.numeric(R.version$minor)
    if (Rmajor == 3 & Rminor == 5) { # we are supporting R >= 3.5.0
        expect_identical(sum(rl$selected), 12L)
    } else {
        expect_identical(sum(rl$selected), 11L)
    }
})

