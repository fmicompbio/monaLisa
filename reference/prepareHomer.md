# Prepare input files for HOMER motif enrichment analysis.

For each bin, write genomic coordinates for foreground and background
regions into files for HOMER motif enrichment analysis.

## Usage

``` r
prepareHomer(
  gr,
  b,
  genomedir,
  outdir,
  motifFile,
  homerfile = findHomer(),
  regionsize = "given",
  Ncpu = 2L,
  verbose = FALSE
)
```

## Arguments

- gr:

  A `GRanges` object (or an object that can be coerced to one) with the
  genomic regions to analyze.

- b:

  A vector of the same length as `gr` that groups its elements into bins
  (typically a factor).

- genomedir:

  Directory containing sequence files in Fasta format (one per
  chromosome).

- outdir:

  A path specifying the folder into which the output files (two files
  per unique value of `b`) will be written.

- motifFile:

  A file with HOMER formatted PWMs to be used in the enrichment
  analysis.

- homerfile:

  Path and file name of the `findMotifsGenome.pl` HOMER script.

- regionsize:

  The peak size to use in HOMER (`"given"` keeps the coordinate region,
  an integer value will keep only that many bases in the region center).

- Ncpu:

  Number of parallel threads that HOMER can use.

- verbose:

  A logical scalar. If `TRUE`, print progress messages.

## Value

The path and name of the script file to run the HOMER motif enrichment
analysis.

## Details

For each bin (unique value of `b`) this functions creates two files in
`outdir` (`outdir/bin_N_foreground.tab` and
`outdir/bin_N_background.tab`, where `N` is the number of the bin and
foreground/background correspond to the ranges that are/are not within
the current bin). The files are in the HOMER peak file format (see
http://homer.ucsd.edu/homer/ngs/peakMotifs.html for details).

In addition, a shell script file is created containing the shell
commands to run the HOMER motif enrichment analysis.

## Examples

``` r
# prepare genome directory (here: one dummy chromosome)
genomedir <- tempfile()
dir.create(genomedir)
writeLines(c(">chr1", "ATGCATGCATCGATCGATCGATCGTACGTA"),
           file.path(genomedir, "chr1.fa"))

# prepare motif file, regions and bins
motiffile <- tempfile()
dumpJaspar(filename = motiffile, pkg = "JASPAR2020",
           opts = list(ID = c("MA0006.1")))
#> [1] TRUE
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:4, width = 4))
b <- bin(1:4, nElements = 2)

# create dummy file (should point to local Homer installation)
homerfile <- file.path(tempdir(), "findMotifsGenome.pl")
writeLines("dummy", homerfile)

# run prepareHomer
outdir <- tempfile()
prepareHomer(gr = gr, b = b, genomedir = genomedir,
             outdir = outdir, motifFile = motiffile,
             homerfile = homerfile, verbose = TRUE)
#> ℹ creating foreground/background region files for HOMER
#> ✔ creating foreground/background region files for HOMER [5ms]
#> 
#> ℹ bin [1,2.5]
#> ✔ bin (2.5,4] [14ms]
#> 
#> ℹ bin (2.5,4]
#> ✔ bin (2.5,4] [17ms]
#> 
#> [1] "/var/folders/8d/778wjbv96mq1760tv6gk374m0000gn/T//RtmpyV8oRz/filee7f739b27b12/run.sh"
list.files(outdir)
#> [1] "bin_001_background.tab" "bin_001_foreground.tab" "bin_002_background.tab"
#> [4] "bin_002_foreground.tab" "run.sh"                

# clean up example
unlink(c(genomedir, motiffile, homerfile, outdir))
```
