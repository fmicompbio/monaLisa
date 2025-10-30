# Adjust for k-mer composition (single iteration)

Adjust background sequence weights for differences in k-mer composition
compared to the foreground sequences. This function implements a single
iteration, and is called iteratively by `.iterativeNormForKmers` to get
to the final set of adjusted weights, which will be the result of
adjusting for GC and k-mer composition. The logic is based on Homer's
`normalizeSequenceIteration()` function found in `Motif2.cpp`.

## Usage

``` r
.normForKmers(
  kmerFreq,
  goodKmers,
  kmerRC,
  seqWgt,
  isForeground,
  minSeqWgt = 0.001,
  maxSeqWgt = 1000
)
```

## Arguments

- kmerFreq:

  A `list` with of matrices. The matrix at index `i` in the list
  contains the probability of k-mers of length `i`, for each k-mer
  (columns) and sequence (rows).

- goodKmers:

  A `list` of `numeric` vectors; the element at index `i` contains the
  number of good (non-N-containing) k-mers of length `i` for each
  sequence.

- kmerRC:

  A `list` of character vectors; the element at index `i` contains the
  reverse complement sequences of all k-mers of length `i`.

- seqWgt:

  A `numeric` vector with starting sequence weights at the beginning of
  the iteration.

- isForeground:

  Logical vector of the same length as `seqs`. `TRUE` indicates that the
  sequence is from the foreground, `FALSE` that it is a background
  sequence.

- minSeqWgt:

  Numeric scalar greater than zero giving the minimal weight of a
  sequence. The default value (0.001) is based on `Homer`
  (HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).

- maxSeqWgt:

  Numeric scalar greater than zero giving the maximal weight of a
  sequence. The default value (1000) is based on `HOMER` (1 /
  HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).

## Value

A named `list` with elements `seqWgt` (updated weights) and `err` (error
measuring difference of foreground and weighted background sequence
compositions).
