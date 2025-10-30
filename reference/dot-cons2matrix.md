# Create matrix from consensus sequence

Given a nucleotide sequence of A,C,G,T letter corresponding to a motif's
consensus string, construct a positional frequency matrix. This matrix
can for example be used as the `profileMatrix` argument in the
constructor for a
[`TFBSTools::PFMatrix`](https://rdrr.io/pkg/TFBSTools/man/XMatrix-class.html)
object.

## Usage

``` r
.cons2matrix(x, n = 100L)
```

## Arguments

- x:

  Character scalar with the motif the consensus sequence.

- n:

  Integer scalar giving the columns sums in the constructed matrix
  (number of observed bases at each position).

## Value

A positional frequency matrix.
