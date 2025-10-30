# Randomized Lasso

This function performs randomized lasso using the `glmnet` package. The
function present in the `stabs` package that runs the lasso version was
adapted for the randomized lasso here. Randomized lasso stability
selection uses this function repeatedly to select predictors.

## Usage

``` r
.glmnetRandomizedLasso(
  x,
  y,
  q,
  weakness = 1,
  type = c("conservative", "anticonservative"),
  ...
)
```

## Arguments

- x:

  The predictor matrix. Passed to `x` of `glmnet.lasso` from `stabs`
  package.

- y:

  The response vector. Passed to `y` of `glmnet.lasso` from `stabs`
  package.

- q:

  The number of variables that are selected on each subsample. Passed to
  `q` of `glmnet.lasso` from `stabs` package.

- weakness:

  Weakness parameter used in randomized lasso (see details).

- type:

  Parameter passed to `type` of `glmnet.lasso` from `stabs` package. It
  is a character vector specifying how much the PFER should be
  controlled. If type is "conservative" (default), then the number of
  selected variables per subsample is \<= q. If type is
  "anticonservative" then the number of selected variables per subsample
  is \>= q.

- ...:

  Additional parameters for `glmnet`.

## Value

The regression output which consists of a list of length 2. The list
contains the following:

- selected:

  \- a logical vector of length equal to the total number of predictors.
  The predictors that were chosen have a value of TRUE.

- path:

  \- a logical matrix containing the regularization steps as columns and
  the predictors as rows. An entry of TRUE indicates selection.

## Details

This function is identical to `glmnet.lasso` from the `stabs` package.
The only addition/modification is the weakness parameter which has been
added when calling the `glmnet` function by setting penalty.factor =
1/runif(ncol(x), weakness, 1), where ncol(x) is the number of
predictors.

## See also

[`glmnet.lasso`](https://rdrr.io/pkg/stabs/man/fitfuns.html) and
[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)
