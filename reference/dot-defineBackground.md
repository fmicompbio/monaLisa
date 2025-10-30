# Define background sequence set for a single motif enrichment calculation

Define the background set for the motif enrichment calculation in a
single bin, depending on the background mode and given foreground
sequences.

## Usage

``` r
.defineBackground(
  sqs,
  bns,
  bg,
  currbn,
  gnm,
  gnm.regions,
  gnm.oversample,
  maxFracN = 0.7,
  GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)
)
```

## Arguments

- sqs, bns, bg:

  The `seqs`, `bins` and `background` arguments from
  `calcBinnedMotifEnrR`.

- currbn:

  An `integer` scalar with the current bin defining the foreground
  sequences.

- gnm, gnm.regions, gnm.oversample:

  The `genome`, `genome.regions` and `genome.oversample` arguments from
  `calcBinnedMotifEnrR`.

- maxFracN:

  The `maxFracN` argument from `calcBinnedMotifEnrR`.

- GCbreaks:

  The breaks between GC bins. The default value is based on the
  hard-coded bins used in Homer.

## Value

A `DataFrame` with sequences represented by rows and columns `seqs`,
`isForeground`, `GCfrac`, `GCbin`, `GCwgt` and `seqWgt`. Only the first
three are already filled in.
