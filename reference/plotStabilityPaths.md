# Plot Stability Paths

Plot the stability paths of each variable (predictor), showing the
selection probability as a function of the regularization step.

## Usage

``` r
plotStabilityPaths(
  se,
  selProbMin = metadata(se)$stabsel.params.cutoff,
  selColor = "cadetblue",
  notSelColor = "grey",
  selProbCutoffColor = "firebrick",
  linewidth = 0.5,
  alpha = 1,
  ylim = c(0, 1),
  labelPaths = FALSE,
  labels = NULL,
  labelNudgeX = 8,
  labelSize = 3
)
```

## Arguments

- se:

  The `SummarizedExperiment` object resulting from stability selection,
  by running
  [`randLassoStabSel`](https://fmicompbio.github.io/monaLisa/reference/randLassoStabSel.md).

- selProbMin:

  A numerical scalar in \[0,1\]. Predictors with a selection probability
  greater than `selProbMin` are shown as colored lines. The color is
  defined by the `col` argument.

- selColor:

  Color for the selected predictors which have a selection probability
  greater than `selProbMin`.

- notSelColor:

  Color for the rest of the (un-selected) predictors.

- selProbCutoffColor:

  Color for the line depicting the selection probability cutoff.

- linewidth:

  Line width.

- alpha:

  Line transparency of the stability paths.

- ylim:

  Limits for y-axis.

- labelPaths:

  If `TRUE`, the predictor labels will be shown at the end of the
  stability paths. The predictor labels given in `labels` will be shown.
  If unspecified, the labels corresponding to the selected predictors
  will be added. If predictors have the same y-value in the last
  regularization step, the labels will be shown in a random order. One
  needs to use `set.seed` to reproduce the plot in this case.

- labels:

  If `labelPaths = TRUE`, the predictors which should be labelled. If
  `NULL`, the selected predictors greater than
  `metadata(se)$stabsel.params.cutoff` will be shown.

- labelNudgeX:

  If `labelPaths = TRUE`, how much to nudge the labels to the right of
  the x-axis.

- labelSize:

  If `labelPaths = TRUE`, the size of the labels.

## Value

a `ggplot2` object.

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
plotStabilityPaths(ss)

```
