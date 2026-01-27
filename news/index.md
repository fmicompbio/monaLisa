# Changelog

## monaLisa 1.17.1

- clarifications of supported families for the stability selection

## monaLisa 1.15.3

- minor updates to adjust to breaking changes in `ggplot2` version 4

## monaLisa 1.15.2

- minor updates to adjust to functionality migrated from `GenomeInfoDb`
  to new `Seqinfo` package

## monaLisa 1.13.3

- switched to `ggplot2` with all plotting functions (except
  `plotMotifHeatmaps` that is still based on `ComplexHeatmaps`)
- new arguments introduced in
  [`plotStabilityPaths()`](https://fmicompbio.github.io/monaLisa/reference/plotStabilityPaths.md)
  to allow for labels to be shown at the end of the stability paths if
  desired
- changed some argument names in
  [`plotStabilityPaths()`](https://fmicompbio.github.io/monaLisa/reference/plotStabilityPaths.md)
  and
  [`plotStabilityPaths()`](https://fmicompbio.github.io/monaLisa/reference/plotStabilityPaths.md).
  See the help pages of these function for more details
- updated the stability selection vignette to showcase labeling
  predictors
- exposed `glmnet` arguments in
  [`randLassoStabSel()`](https://fmicompbio.github.io/monaLisa/reference/randLassoStabSel.md)

## monaLisa 1.13.2

- minor fix to ensure that heatmap bin legends are ordered in the same
  way as the heatmap columns

## monaLisa 1.11.3

- add show_bin_legend argument to plotMotifHeatmaps (contributed by
  [@danymukesha](https://github.com/danymukesha), PR
  [\#62](https://github.com/fmicompbio/monaLisa/issues/62))

## monaLisa 1.9.1

- adapt `dumpJaspar` to also work with `JASPAR2024`

## monaLisa 1.7.1

- allow modification of heatmap graphical parameters by forwarding `...`
  argument in `plotMotifHeatmaps` to all calls to
  [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## monaLisa 1.3.1

- update citation information

## monaLisa 1.1.2

- added citation to the Bioinformatics publication

## monaLisa 1.1.1

- added link to pre-print manuscript on biorXiv to `README.md`
- added warning to `bin(..., minAbsX = val)` if adjusted zero-bin breaks
  deviate more than 20% from `val`
- added `doPlot` argument to `plotMotifHeatmaps` to select if heatmaps
  should be plotted or just generated and returned
- added `LICENSE.md` file
- expanded `monaLisa.Rmd` vignette with illustration on how to do a
  binary or single set motif enrichment analysis
- expanded on collineairty in regression in the
  `selecting_motifs_with_randLassoStabSel.Rmd` vignette and the choice
  of parameter values in stability selection.
- updated the `results.binned_6mer_enrichment_LMRs.rds` and
  `results.binned_motif_enrichment_LMRs.rds` files stored in `monaLisa`
  under the current version of the package.

## monaLisa 1.0.0

- Initial release of `monaLisa` as part of Bioconductor 3.14

## monaLisa 0.99.5

- Updated `R/monaLisa-package.R` file

## monaLisa 0.99.4

- Suppressed warnings from matchPWM (due to presence of Ns) in
  regression vignette

## monaLisa 0.99.3

- Updated `README.md` file

## monaLisa 0.99.2

- Added fixes to the regression vignette
- Addressed failing test in `calcBinnedKmerEnr`

## monaLisa 0.99.1

- Added examples where missing for exported functions
- Harmonized function naming (anno_seqlogo -\> annoSeqlogo,
  sample_random_regions -\> sampleRandomRegions)
- Clarified details on Pearson residual calculation
- Adapted documentation for new version of BiocParallel
- Harmonized return values from plot functions
- Added legend position and size arguments to plotSelectionProb()

## monaLisa 0.99.0

- Preparation for Bioconductor submission

## monaLisa 0.2.0

- Added a `NEWS.md` file to track changes to the package.
