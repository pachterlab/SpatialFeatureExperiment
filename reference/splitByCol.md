# Split SFE object with categorical vector or geometry

The `split` methods for SFE split an SFE object into multiple SFE
objects by geometries (all cells/spots intersecting with each geometry
will become a separate SFE object). The `splitSamples` function splits
the SFE object by `sample_id` so each sample will become a separate SFE
object. The `splitContiguity` function splits the SFE object by
contiguity of an `annotGeometry`, which by default is "tissueBoundary".
The `splitComponent` function splits the SFE object by graph component,
so if there are disconnected components in the graph, then each
component will become a new SFE object.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment,sf'
splitByCol(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE)

# S4 method for class 'SpatialFeatureExperiment,sfc'
splitByCol(x, f, sample_id = 1L, colGeometryName = 1L, cover = FALSE)

# S4 method for class 'SpatialFeatureExperiment,list'
splitByCol(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE)

splitSamples(x)

splitContiguity(
  x,
  colGeometryName = 1L,
  annotGeometryName = "tissueBoundary",
  min_area = 0,
  cover = FALSE
)
```

## Arguments

- x:

  An SFE object

- f:

  It can be a `sf` data frame or `sfc` to split by geometry. Each row of
  the `sf` data frame or each element in the `sfc` will correspond to a
  new SFE object. The `sf` data frame must have a column `sample_id`
  when splitting multiple samples. Can also be a list of `sfc` whose
  names correspond to `sample_id`s to split.

- sample_id:

  Which samples to split.

- colGeometryName:

  Which `colGeometry` to use to determine which cells or spots should
  belong to which new SFE object when splitting by `sf` or `sfc`.
  Default to the first one.

- cover:

  Logical, whether the geometries in `x` must be entirely covered by `y`
  if `op = st_intersection` or whether `x` must be entirely outside `y`
  if `op = st_difference`. Only relevant when `keep_whole != "none"`.

- annotGeometryName:

  Name of `annotGeometry` to use to split by contiguity.

- min_area:

  Minimum area in the same unit as the geometry coordinates (squared)
  for each piece to be considered a separate piece when splitting by
  contiguity. Only pieces that are large enough are considered.

## Value

A list of SFE objects.

## Examples

``` r
# example code
```
