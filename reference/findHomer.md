# Find HOMER script file.

Find absolute path to HOMER script file.

## Usage

``` r
findHomer(homerfile = "findMotifsGenome.pl", dirs = NULL)
```

## Arguments

- homerfile:

  Name of the script file to search.

- dirs:

  Directory names to look for `homerfile`. If `dirs=NULL`, all
  directories listed in the `PATH` environment variable will be
  searched.

## Value

Absolute path to `homerfile`, or `NA` if none or several were found.

## Details

In addition to `dirs`, `findHomer` will also look in the directory
provided in the environment variable `MONALISA_HOMER`.

## Examples

``` r
homer_path <- findHomer()
```
