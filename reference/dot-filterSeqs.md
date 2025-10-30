# Filter Sequences

Filter sequences that are unlikely to be useful for motif enrichment
analysis. The current defaults are based on HOMER (version 4.11).

## Usage

``` r
.filterSeqs(
  seqs,
  maxFracN = 0.7,
  minLength = 5L,
  maxLength = 100000L,
  verbose = FALSE
)
```

## Arguments

- seqs:

  A `DNAStringSet` object.

- maxFracN:

  A numeric scalar with the maximal fraction of N bases allowed in a
  sequence (defaults to 0.7).

- minLength:

  The minimum sequence length (default from Homer). Sequences shorter
  than this will be filtered out.

- maxLength:

  The maximum sequence length (default from Homer). Sequences bigger
  than this will be filtered out.

- verbose:

  A logical scalar. If `TRUE`, report on filtering.

## Value

A logical vector of the same length as `seqs` with `TRUE` indicated to
keep the sequence and `FALSE` to filter it out.

## Details

The filtering logic is based on `removePoorSeq.pl` from Homer.
