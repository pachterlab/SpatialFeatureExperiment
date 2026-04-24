# Binary predicates for geometry of each cell/spot and annotation

This function finds binary predicates for the geometry of each cell/spot
(i.e. `colGeometry`) and an annotation geometry for each sample. For
example, whether each Visium spot intersects with the tissue boundary in
each sample.

## Usage

``` r
annotPred(
  sfe,
  colGeometryName = 1L,
  annotGeometryName = 1L,
  sample_id = "all",
  pred = st_intersects,
  yx = FALSE
)

annotNPred(
  sfe,
  colGeometryName = 1L,
  annotGeometryName = 1L,
  sample_id = "all",
  pred = st_intersects
)
```

## Arguments

- sfe:

  An SFE object.

- colGeometryName:

  Name of column geometry for the predicate.

- annotGeometryName:

  Name of annotation geometry for the predicate.

- sample_id:

  Which sample(s) to operate on. Can be "all" to indicate all samples.

- pred:

  Predicate function to use, defaults to
  [`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

- yx:

  Whether to do `pred(y, x)` instead of `pred(x, y)`. For symmetric
  predicates, the results should be the same. When x has a large number
  of geometries and y has few, `pred(y, x)` is much faster than
  `pred(x, y)` for `st_intersects`, `st_disjoint`, and
  `st_is_within_distance`.

## Value

For `annotPred`, a logical vector of the same length as the number of
columns in the sample(s) of interest, with barcodes (or corresponding
column names of sfe) as names. For `annotNPred`, a numeric vector of the
same length as the number of columns in the sample(s) of interest with
barcodes as names, indicating the number of geometries in the
`annotGeometry` of interest returns TRUE for the predicate for each each
geometry in the `colGeometry` of interest.

## See also

annotOp

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
# Whether each spot is in tissue
in_tissue <- annotPred(sfe, "spotPoly", annotGeometryName = "tissueBoundary")
# How many nuclei are there in each Visium spot
n_nuclei <- annotNPred(sfe, "spotPoly", annotGeometryName = "nuclei")
```
