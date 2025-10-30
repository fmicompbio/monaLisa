# Check if elements of \`x\` are have equal lengths

Check if the elements of \`x\` are all equally long. If not, generate a
warning.

## Usage

``` r
.checkIfSeqsAreEqualLength(x)
```

## Arguments

- x:

  An object that implements a `width` method, typically a `GRanges` or
  `DNAStringSet` object.

## Value

`NULL` (invisibly). The function is called for its side-effect of
generating a warning if elements of the input are not of equal lengths.
