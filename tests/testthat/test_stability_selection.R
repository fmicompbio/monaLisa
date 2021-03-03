test_that("randLassoStabSel() works properly", {

  # create data set
  set.seed(555)
  Y <- rnorm(n = 500, mean = 2, sd = 1)
  X <- matrix(data = NA, nrow = length(Y), ncol = 50)
  for (i in 1:ncol(X)) {
    X[ ,i] <- runif(n = 500, min = 0, max = 3)
  }
  s_cols <- sample(x = 1:ncol(X), size = 10, replace = FALSE)
  for (i in 1:length(s_cols)) {
    X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
  }

  X2 <- X
  Y2 <- Y
  rownames(X2) <- paste0("peak", 1:nrow(X2))
  colnames(X2) <- paste0("motif", 1:ncol(X2))
  names(Y2) <- paste0("peak", 1:length(Y2))

  # randomized lasso stability selection
  ss <- monaLisa::randLassoStabSel(x = X, y = Y)

  # tests
  expect_true(is(ss, "SummarizedExperiment"))
  expect_true(all(rowData(ss)$y == Y))
  expect_identical(ss$selProb, colData(ss)[, ncol(colData(ss))])
  expect_true(all(s_cols %in% metadata(ss)$stabsel.params.selected))
  expect_identical(dim(ss), c(500L, 50L))
  expect_identical(length(Y), nrow(ss))
  expect_true(!is.null(SummarizedExperiment::assay(ss)))
  expect_error(randLassoStabSel(x = as.data.frame(X), y = Y))
  expect_error(randLassoStabSel(x = X2[1:100, ], y = Y2[2:101]))
  expect_error(randLassoStabSel(x = X2[1:100, ], y = Y2[2:100]))
})


test_that(".glmnetRandomizedLasso() works properly", {

  # create data set
  set.seed(555)
  Y <- rnorm(n = 500, mean = 2, sd = 1)
  X <- matrix(data = NA, nrow = length(Y), ncol = 50)
  for (i in 1:ncol(X)) {
    X[ ,i] <- runif(n = 500, min = 0, max = 3)
  }
  s_cols <- sample(x = 1:ncol(X), size = 10, replace = FALSE)
  for (i in 1:length(s_cols)) {
    X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
  }

  # tests

  # ... x as data.frame
  expect_message(.glmnetRandomizedLasso(x = as.data.frame(X), y = Y, q=11), "coerced to a model matrix without intercept")

  # ... with type="anticonservative"
  rl <- .glmnetRandomizedLasso(x = X, y = Y, q=11, type="anticonservative")

  # ... with specific lambda
  expect_error(.glmnetRandomizedLasso(x = X, y = Y, q=11, lambda=5))

  # ... expect_true
  expect_true(base::inherits(rl, "list"))
  expect_true(all(names(rl)==c("selected", "path")))
  expect_true(base::inherits(rl$selected, "logical"))
  expect_true(base::inherits(rl$path, "matrix"))

  # outputs differ depending on R version due to random number generator with 'sample' function (check RNGkind) --> different as of R 3.6.0
  Rmajor <- as.numeric(R.version$major)
  Rminor <- as.numeric(R.version$minor)
  if (Rmajor == 3 & Rminor == 5) { # we are supporting R >= 3.5.0
    expect_true(sum(rl$selected)==12)
  } else {
    expect_true(sum(rl$selected)==11)
  }


})

