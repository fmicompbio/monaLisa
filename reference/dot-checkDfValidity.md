# Check if seqinfo DataFrame is valid

Check if the `DataFrame` with sequence information is valid, i.e. is of
the correct object type (`DataFrame`) and has all expected columns and
attributes.

## Usage

``` r
.checkDfValidity(df)
```

## Arguments

- df:

  Input object to be checked. It should have an attribute `err` and
  columns:

  `seqs`

  :   : a `DNAStringSet` object.

  `isForeground`

  :   that indicates if a sequence is in the foreground group.

  `GCfrac`

  :   : the fraction of G+C bases per sequence.

  `GCbin`

  :   : the GC bin for each sequence.

  `GCwgt`

  :   : the sequence weight to adjust for GC differences between
      foreground and background sequences.

  `seqWgt`

  :   : the sequence weight to adjust for k-mer differences between
      foreground and background sequences.

## Value

`TRUE` (invisibly) if `df` is valid, otherwise it raises an exception
using
[`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
