# monaLisa 1.9.1

* adapt `dumpJaspar` to also work with `JASPAR2024`

# monaLisa 1.7.1

* allow modification of heatmap graphical parameters by forwarding `...` argument in `plotMotifHeatmaps` to all calls to `ComplexHeatmap::Heatmap`

# monaLisa 1.3.1

* update citation information

# monaLisa 1.1.2

* added citation to the Bioinformatics publication

# monaLisa 1.1.1

* added link to pre-print manuscript on biorXiv to `README.md`
* added warning to `bin(..., minAbsX = val)` if adjusted zero-bin breaks deviate more than 20% from `val`
* added `doPlot` argument to `plotMotifHeatmaps` to select if heatmaps should be plotted or just generated and returned
* added `LICENSE.md` file
* expanded `monaLisa.Rmd` vignette with illustration on how to do a binary or single set motif enrichment analysis
* expanded on collineairty in regression in the `selecting_motifs_with_randLassoStabSel.Rmd` vignette and the choice of parameter values in stability selection.
* updated the `results.binned_6mer_enrichment_LMRs.rds` and `results.binned_motif_enrichment_LMRs.rds` files stored in `monaLisa` under the current version of the package.

# monaLisa 1.0.0

* Initial release of `monaLisa` as part of Bioconductor 3.14

# monaLisa 0.99.5

* Updated `R/monaLisa-package.R` file

# monaLisa 0.99.4

* Suppressed warnings from matchPWM (due to presence of Ns) in regression vignette

# monaLisa 0.99.3

* Updated `README.md` file

# monaLisa 0.99.2

* Added fixes to the regression vignette
* Addressed failing test in `calcBinnedKmerEnr`

# monaLisa 0.99.1

* Added examples where missing for exported functions
* Harmonized function naming (anno_seqlogo -> annoSeqlogo, sample_random_regions -> sampleRandomRegions)
* Clarified details on Pearson residual calculation
* Adapted documentation for new version of BiocParallel
* Harmonized return values from plot functions
* Added legend position and size arguments to plotSelectionProb()

# monaLisa 0.99.0

* Preparation for Bioconductor submission

# monaLisa 0.2.0

* Added a `NEWS.md` file to track changes to the package.
