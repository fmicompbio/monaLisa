# Sample random regions of fixed length.

Sample random regions from the mappable parts of the genome with a given
fraction from CpG islands.

## Usage

``` r
sampleRandomRegions(allowedRegions = NULL, N = 100L, regWidth = 200L)
```

## Arguments

- allowedRegions:

  An unstranded `GRanges` object of the "allowed" of the genome, usually
  the mappable regions.

- N:

  Number of regions to sample.

- regWidth:

  Region width.

## Value

A `GRanges` object with randomly sampled mappable regions of width
`regWidth` with `fractionCGI` coming from CpG islands.

## Details

In order to make the results deterministic, set the random number seed
before calling `sampleRandomRegions` using `set.seed`.

## Examples

``` r
regs <- GenomicRanges::GRanges(
  seqnames = rep(c("chr1", "chr2"), each = 2), 
  ranges = IRanges::IRanges(start = 1:4, end = 5:8))
set.seed(123)
sampleRandomRegions(regs, N = 2, regWidth = 3L)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       3-5      *
#>   [2]     chr2       4-6      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
