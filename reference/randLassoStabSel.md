# Randomized Lasso Stability Selection

This function runs randomized lasso stability selection as presented by
Meinshausen and Bühlmann (2010) and with the improved error bounds
introduced by Shah and Samworth (2013). The function uses the
[`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html) function from
the `stabs` package, but implements the randomized lasso version.

## Usage

``` r
randLassoStabSel(
  x,
  y,
  weakness = 0.8,
  cutoff = 0.8,
  PFER = 2,
  mc.cores = 1L,
  glmnet.args = list(),
  ...
)
```

## Arguments

- x:

  The predictor matrix.

- y:

  The response vector. This should be a numeric vector (also for the
  binomial case - in this case it will be converted to a factor
  internally by `glmnet`).

- weakness:

  Value between 0 and 1 (default = 0.8). It affects how strict the
  method will be in selecting predictors. The closer it is to 0, the
  more stringent the selection. A weakness value of 1 is identical to
  performing lasso stability selection (not the randomized version).

- cutoff:

  Value between 0 and 1 (default = 0.8) which is the cutoff for the
  selection probability. Any variable with a selection probability that
  is higher than the set cutoff will be selected.

- PFER:

  Integer (default = 2) representing the absolute number of false
  positives that we allow for in the final list of selected variables.
  For details see Meinshausen and Bühlmann (2010).

- mc.cores:

  Integer (default = 1) specifying the number of cores to use in
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html), which is the
  default way [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)
  does parallelization.

- glmnet.args:

  Named list with additional arguments to the internal
  `.glmnetRandomizedLasso` function (beyond `x`, `y` and `weakness`,
  which are determined automatically, and `q`, which should not be
  specified (it will be determined from `cutoff` and `PFER`). The
  available arguments to `.glmnetRandomizedLasso` are the same as the
  ones for [`glmnet.lasso`](https://rdrr.io/pkg/stabs/man/fitfuns.html).
  A typical use case would be to define the `family` argument to
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)
  (currently "gaussian" and "binomial" are supported).

- ...:

  Additional parameters that can be passed on to
  [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html).

## Value

A `SummarizedExperiment` object where the rows are the observations and
the columns the predictors (same dimnames as the predictor matrix `x`).
It contains:

- assays:

  :

  x

  :   : the predictor matrix.

- rowData:

  : a `DataFrame` with columns:

  y

  :   : the response vector.

- colData:

  : a `DataFrame` with columns:

  selProb

  :   : the final selection probabilities for the predictors (from the
      last regularization step).

  selected

  :   : logical indicating the predictors that made the selection with
      the specified cutoff.

  selAUC

  :   : the normalized area under the seletion curve (mean of selection
      probabilities over regulatization steps).

  reg'`i`'

  :   : columns containing the selection probabilities for
      regularization step i.

- metadata:

  : a list of output returned from
  [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html) and
  `randLassoStabSel`:

  stabsel.params.cutoff

  :   : probability cutoff set for selection of predictors (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.selected

  :   : elements with maximal selection probability greater `cutoff`
      (see [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.max

  :   : maximum of selection probabilities (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.q

  :   : average number of selected variables used (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.PFER

  :   : (realized) upper bound for the per-family error rate (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.specifiedPFER

  :   : specified upper bound for the per-family error rate (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.p

  :   : the number of effects subject to selection (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.B

  :   : the number of subsamples (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.sampling.type

  :   : the sampling type used for stability selection (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.assumption

  :   : the assumptions made on the selection probabilities (see
      [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)).

  stabsel.params.call

  :   : [`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html) the
      call.

  randStabsel.params.weakness

  :   : the weakness parameter in the randomized lasso stability
      selection.

## Details

Randomized lasso stability selection runs a randomized lasso regression
several times on subsamples of the response variable and predictor
matrix. N/2 elements from the response variable are randomly chosen in
each regression, where N is the length of the vector. The corresponding
section of the predictor matrix is also chosen, and the internal
`.glmnetRandomizedLasso` function is applied. Stability selection
results in selection probabilities for each predictor. The probability
of a specific predictor is the number of times it was selected divided
by the total number of subsamples that were done (total number of times
the regression was performed).

We made use of the `stabs` package that implements lasso stability
selection, and adapted it to run randomized lasso stability selection.

## References

N. Meinshausen and P. Bühlmann (2010), Stability Selection, *Journal of
the Royal Statistical Society: Series B (Statistical Methodology)*,
**72**, 417–73.  
R.D. Shah and R.J. Samworth (2013), Variable Selection with Error
Control: Another Look at Stability Selection, *Journal of the Royal
Statistical Society: Series B (Statistical Methodology)*, **75**,
55–80.  
B. Hofner, L. Boccuto, and M. Göker (2015), Controlling False
Discoveries in High-Dimensional Situations: Boosting with Stability
Selection, *BMC Bioinformatics*, **16** 144.

## See also

[`stabsel`](https://rdrr.io/pkg/stabs/man/stabsel.html)

## Examples

``` r
## create data set
Y <- rnorm(n = 500, mean = 2, sd = 1)
X <- matrix(data = NA, nrow = length(Y), ncol = 50)
for (i in seq_len(ncol(X))) {
  X[ ,i] <- runif(n = 500, min = 0, max = 3)
}
s_cols <- sample(x = seq_len(ncol(X)), size = 10,
  replace = FALSE)
for (i in seq_along(s_cols)) {
  X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
}

## reproducible randLassoStabSel() with 1 core
set.seed(123)
ss <- randLassoStabSel(x = X, y = Y)

## reproducible randLassoStabSel() in parallel mode
## (only works on non-windows machines)
if (.Platform$OS.type == "unix") {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
    ss <- randLassoStabSel(x = X, y = Y, mc.preschedule = TRUE,
                           mc.set.seed = TRUE, mc.cores = 2L)
}
```
