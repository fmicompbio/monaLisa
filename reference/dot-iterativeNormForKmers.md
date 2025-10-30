# Adjust for k-mer composition (multiple iterations)

Here we run \`.normForKmers\` multiple times to converge to the final
weights that will be used to correct the background sequences for k-mer
composition differences compared to the foreground. We closely follow
`HOMER`'s `normalizeSequence()` function found in `Motif2.cpp`. Note
that `HOMER` runs the `normalizeSequence()` one last time after going
through all iterations or reaching a low error, which we do not do here.

## Usage

``` r
.iterativeNormForKmers(
  df,
  maxKmerSize = 3L,
  minSeqWgt = 0.001,
  maxIter = 160L,
  verbose = FALSE
)
```

## Arguments

- df:

  A `DataFrame` with sequence information as returned by
  `.calculateGCweight`.

- maxKmerSize:

  Integer scalar giving the maximum k-mer size to consider. The default
  is set to 3 (like in `HOMER`), meaning that k-mers of size 1, 2 and 3
  are considered.

- minSeqWgt:

  Numeric scalar greater than zero giving the minimal weight of a
  sequence. The default value (0.001) was also used by `HOMER`
  (HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).

- maxIter:

  An integer scalar giving the maximum number if times to run
  `.normForKmers`. the default is set to 160 (as in `HOMER`).

- verbose:

  A logical scalar. If `TRUE`, report on k-mer composition adjustment.

## Value

A `DataFrame` containing:

- sequenceWeights:

  : a `dataframe` containing the sequence GC content, GC bins they were
  assigned to, the weight to correct for GC differences between
  foreGround and background sequences, the weight to adjust for kmer
  composition, and the the error term

- sequenceNucleotides:

  : a `DNAStringSet` object containing the raw sequences
