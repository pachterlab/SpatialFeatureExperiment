# Binary operations for geometry of each cell/spot and annotation

**\[deprecated\]** Just like
[`annotPred`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotPred.md),
but performs the operation rather than predicate. For example, this
function would return the geometry of the intersections between each
Visium spot and the tissue boundary for each sample, rather than whether
each Visium spot intersects the tissue boundary. In case one cell/spot
gets broken up into multiple geometries, the union of those geometries
will be taken, so each cell/spot will only get one geometry.

## Usage

``` r
annotOp(
  sfe,
  colGeometryName = 1L,
  annotGeometryName = 1L,
  sample_id = "all",
  op = deprecated(),
  difference = FALSE
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

- op:

  Deprecated.

- difference:

  Logical, whether to compute remove the area occupied by `y` from the
  SFE object.

## Value

A `sf` data frame with `geometry` column containing the geometries and
corresponding column names of sfe as row names. There is no guarantee
that the returned geometries are valid or preserve the geometry class
(e.g. when the intersection of polygons result into a line of a point).

## See also

annotPred

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
# Get the intersection of myofibers with each Visium spot
myofibers_on_spots <- annotOp(sfe, "spotPoly",
    annotGeometryName = "myofiber_simplified"
)
#> Warning: `annotOp()` was deprecated in SpatialFeatureExperiment 1.16.0.
#> ℹ In addition, DuckDBDataFrame is not supported
```
