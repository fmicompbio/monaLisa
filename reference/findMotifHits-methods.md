# Find motif matches in sequences.

`findMotifHits` scans sequences (either provided as a file, an R object
or genomic coordinates) for matches to positional weight matrices
(provided as a file or as R objects)

## Usage

``` r
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'character,character'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'character,DNAString'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'character,DNAStringSet'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrix,character'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrix,DNAString'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrix,DNAStringSet'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrixList,character'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrixList,DNAString'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrixList,DNAStringSet'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrix,GRanges'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)

# S4 method for class 'PWMatrixList,GRanges'
findMotifHits(
  query,
  subject,
  min.score,
  method = c("matchPWM", "homer2"),
  homerfile = findHomer("homer2"),
  BPPARAM = SerialParam(),
  genome = NULL
)
```

## Arguments

- query:

  The motifs to search for, either a

  `character(1)`

  :   with the path and file name of a motif file with PWM in HOMER
      format (currently only supported for `method="homer2"`)

  `PWMatrix`

  :   with a single PWM

  `PWMatrixList`

  :   with several PWMs to search for.

- subject:

  The sequences to be searched, either a

  `character`

  :   with the path and file name of a sequence file with DNA sequences
      in FASTA format

  `DNAString`

  :   with a single sequence

  `DNAStringSet`

  :   with several sequences

  `GRanges`

  :   object with the genomic coordinates of the sequences to be
      searched.

- min.score:

  The minimum score for counting a match. Can be given as a character
  string containing a percentage (e.g. "85 highest possible score or as
  a single number.

- method:

  The internal method to use for motif searching. One of

  `"matchPWM"`

  :   using Biostrings::matchPWM (optimized)

  `"homer2"`

  :   call to the homer2 binary

  Please note that the two methods might give slightly different results
  (see details).

- homerfile:

  Path and file name of the `homer2` binary.

- BPPARAM:

  An optional
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance determining the parallel back-end to be used during
  evaluation.

- genome:

  `BSgenome` object that is the reference genome of the subject. This
  argument is set to NULL by default and only used by the function when
  the subject is a `GRanges` object. It is then necessary to specify the
  genome so that the function can internally convert the genomic regions
  into a `DNAStringSet` object.

## Value

A `GRanges` object with the matches to `query` in `subject`.

## Details

The implemented methods (`matchPWM` and `homer2`) are there for
convenience (`method="matchPWM"` calls
[`Biostrings::matchPWM`](https://rdrr.io/pkg/Biostrings/man/matchPWM.html)
internally in an optimized fashion, and `method = "homer2"` calls the
command line tool from Homer and therefore requires an installation of
Homer).

In general, running `findMotifHits` with the same parameters using any
of the methods generates identical results. Some minor differences could
occur that result from rounding errors during the necessary conversion
of PWMs (log2-odd scores) to the probability matrices needed by Homer,
and the conversion of scores from and to the natural log scale used by
Homer. These conversions are implemented transparently for the user, so
that the arguments of `findMotifHits` do not have to be adjusted (e.g.
the PWMs should always contain log2-odd scores, and `min.score` is
always on the log2 scale).

If there are bases with frequencies of less than 0.001 in a motif, Homer
will set them to 0.001 and adjust the other frequencies at that motif
position accordingly so that they sum to 1.0. This may differ from the
adjustment used when scanning a PWM with `matchPWM` (e.g. the
`pseudocounts` argument in the
[`toPWM`](https://rdrr.io/pkg/TFBSTools/man/toPWM-methods.html)
function), and thus can give rise to differences in reported motif hits
and hit scores (typically only low-scoring hits).

## Examples

``` r
seqs <- Biostrings::DNAStringSet(c(s1 = "GTCAGTCGATC", s2 = "CAGTCTAGCTG",
                                   s3 = "CGATCGTCAGT", s4 = "AGCTGCAGTCT"))
m <- rbind(A = c(2, 0, 0),
           C = c(1, 1, 0),
           G = c(0, 2, 0),
           T = c(0, 0, 3))
pwms <- TFBSTools::PWMatrixList(
    TFBSTools::PWMatrix(ID = "m1", profileMatrix = m),
    TFBSTools::PWMatrix(ID = "m2", profileMatrix = m[, 3:1])
)
findMotifHits(pwms, seqs, min.score = 7)
#> GRanges object with 6 ranges and 4 metadata columns:
#>       seqnames    ranges strand |     matchedSeq pwmid pwmname     score
#>          <Rle> <IRanges>  <Rle> | <DNAStringSet> <Rle>   <Rle> <numeric>
#>   [1]       s1       4-6      + |            AGT    m1 Unknown         7
#>   [2]       s1       2-4      - |            TGA    m2 Unknown         7
#>   [3]       s2       2-4      + |            AGT    m1 Unknown         7
#>   [4]       s3      9-11      + |            AGT    m1 Unknown         7
#>   [5]       s3       7-9      - |            TGA    m2 Unknown         7
#>   [6]       s4       7-9      + |            AGT    m1 Unknown         7
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome
```
