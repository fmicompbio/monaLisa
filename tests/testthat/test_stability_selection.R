context("stability selection")

test_that("randomized_stabsel() works properly", {

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

  # randomized lasso stability selection
  ss <- lisa::randomized_stabsel(x = X, y = Y)

  # tests
  expect_true(base::inherits(ss, "stabsel"))
  expect_true(all(s_cols %in% ss$selected))
  expect_true(all(!(seq(1,ncol(X),1)[-s_cols] %in% ss$selected)))
  expect_true(all(c("phat", "selected", "cutoff", "PFER") %in% names(ss))) # this is crucial because they are accessed in other functions in lisa
  expect_true(!is.null(ss$phat))

})


test_that("glmnet.randomized_lasso() works properly", {
  
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
  expect_message(lisa::glmnet.randomized_lasso(x = as.data.frame(X), y = Y, q=11), "coerced to a model matrix without intercept")
  
  # ... with type="anticonservative"
  rl <- lisa::glmnet.randomized_lasso(x = X, y = Y, q=11, type="anticonservative")
  
  # ... with specific lambda
  expect_error(lisa::glmnet.randomized_lasso(x = X, y = Y, q=11, lambda=5))
  
  # ... expect_true
  expect_true(base::inherits(rl, "list"))
  expect_true(all(names(rl)==c("selected", "path")))
  expect_true(base::inherits(rl$selected, "logical"))
  expect_true(base::inherits(rl$path, "matrix"))
  
  # outputs differ depending on R version due to random number generator with 'sample' function (check RNGkind) --> different as of R 3.6.0
  Rmajor <- as.numeric(R.version$major)
  Rminor <- as.numeric(R.version$minor)
  if (Rmajor>=3 & Rminor>=6) {
    expect_true(sum(rl$selected)==11)
  } else {
    expect_true(sum(rl$selected)==12)
  }

  
})

