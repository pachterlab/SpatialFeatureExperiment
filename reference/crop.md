# Crop an SFE object with a geometry

Returns an SFE object whose specified `colGeometry` returns `TRUE` with
a geometric predicate function (usually intersects) with another
geometry of interest. This can be used to subset an SFE object with a
tissue boundary or histological region polygon, or crop away empty
spaces. After cropping, not only will the cells/spots be subsetted, but
also all geometries will be cropped.

## Usage

``` r
crop(
  x,
  y = NULL,
  colGeometryName = 1L,
  sample_id = "all",
  op = st_intersection,
  keep_whole = "none",
  cover = FALSE
)
```

## Arguments

- x:

  An SFE object.

- y:

  An object of class `sf`, `sfg`, `sfc` with which to crop the SFE
  object, or a bounding box with the format of the output of
  [`bbox,SpatialFeatureExperiment-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/bbox-SpatialFeatureExperiment-method.md).

- colGeometryName:

  Column geometry to used to indicate which cells/spots to keep.

- sample_id:

  Samples to crop. Optional when only one sample is present. Can be
  multiple samples, or "all", which means all samples. For multiple
  samples, `sf` data frame `y` may have column `sample_id` indicating
  which geometry subsets which sample or matrix `y` may indicate sample
  specific bounding boxes in its column names. Only samples included in
  the indicated sample IDs are subsetted. If sample is not indicated in
  `y`, then the same geometry or bounding box is used to subset all
  samples specified in the `sample_id` argument.

- op:

  A geometric operation function to crop the geometries in the SFE
  object. Only
  [`st_intersection`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
  and
  [`st_difference`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
  are allowed. If "intersection", then only things inside `y` is kept
  after cropping. If "difference", then only things outside `y` is kept.

- keep_whole:

  Character vector, can be one or more of "col" and "annot" to keep
  whole items from `colGeometries` or `annotGeometries`, keeping
  geometries that partially intersect with `y` whole. This can greatly
  speed up code while not breaking geometries into multiple pieces. Can
  also be "none" so all geometries are actually cropped.

- cover:

  Logical, whether the geometries in `x` must be entirely covered by `y`
  if `op = st_intersection` or whether `x` must be entirely outside `y`
  if `op = st_difference`. Only relevant when `keep_whole != "none"`.

## Value

An SFE object. There is no guarantee that the geometries after cropping
are still all valid or preserve the original geometry class.

## Details

3D geometries are allowed, but geometric operations can only be
performed in x and y but not z.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
# Subset sfe to only keep spots on tissue
sfe_on_tissue <- crop(sfe, tissueBoundary(sfe),
    colGeometryName = "spotPoly",
    sample_id = "Vis5A"
)
```
