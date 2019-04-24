test_that("randomized_stabsel() works properly", {
  
  # create data set
  set.seed(555)
  Y <- rnorm(n = 500, mean = 2, sd = 1)
  X <- matrix(data = NA, nrow = length(Y), ncol = 50)
  for (i in 1:ncol(X)){
    X[ ,i] <- runif(n = 500, min = 0, max = 3) 
  }
  s_cols <- sample(x = 1:ncol(X), size = 10, replace = FALSE)
  for (i in 1:length(s_cols)){
    X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
  }
  
  # randomized lasso stability selection
  ss <- randomized_stabsel(x=X, y=Y)
  
  # tests
  expect_true(inherits(ss, "stabsel"))
  expect_true(all(s_cols%in%ss$selected))
  expect_true(all(c("phat", "selected", "cutoff", "PFER")%in%names(ss))) # this is crucial because they are accessed in other functions in lisa
  expect_true(!is.null(ss$phat))

})