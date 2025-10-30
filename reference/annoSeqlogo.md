# Sequence logo annotation

Create an annotation for a
[`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
containing sequence logos.

## Usage

``` r
annoSeqlogo(
  grobL,
  which = c("column", "row"),
  space = unit(0.5, "mm"),
  width = NULL,
  height = NULL,
  gp = gpar(fill = NA, col = NA)
)
```

## Arguments

- grobL:

  A `list` of sequence logo grobs, typically created using
  [`seqLogoGrob`](https://fmicompbio.github.io/monaLisa/reference/seqLogoGrob.md).

- which:

  Whether it is a column annotation or a row annotation?

- space:

  The space around the image to the annotation grid borders. The value
  should be a unit object.

- width:

  Width of the annotation. The value should be an absolute unit. Width
  is not allowed to be set for column annotation.

- height:

  Height of the annotation. The value should be an absolute unit. Height
  is not allowed to be set for row annotation.

- gp:

  Graphic parameters for annotation grids. Can be used to control the
  background color in the annotation grids.

## Value

An annotation function which can be used in
[`HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html).

## Examples

``` r
if (require(JASPAR2020) && require(TFBSTools) && require(gridExtra)) {
    pfm1 <- getMatrixByID(JASPAR2020, "MA0139")

    g1 <- seqLogoGrob(pfm1)

    anno <- annoSeqlogo(list(g1))
}
#> Loading required package: JASPAR2020
#> Loading required package: TFBSTools
#> Loading required package: gridExtra
```
