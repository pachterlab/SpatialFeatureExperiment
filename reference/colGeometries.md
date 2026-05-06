# Column geometry getters and setters

`colGeometries` are geometries that correspond to columns of the gene
count matrix, such as Visium spots or cells. Same as
`dimGeometry(x, MARGIN = 2L, ...)`, with convenience wrappers for
getters and setters of special geometries:

- spotPoly:

  Polygons of spots from technologies such as Visium, ST, and slide-seq,
  which do not correspond to cells. Centroids of the polygons are stored
  in `spatialCoords` of the underlying `SpatialExperiment` object.

- ROIPoly:

  Polygons of regions of interest (ROIs) from technologies such as laser
  capture microdissection (LCM) and GeoMX DSP. These should correspond
  to columns of the gene count matrix.

- cellSeg:

  Cell segmentation polygons. If the columns of the gene count matrix
  are single cells, then this is stored in `colGeometries`. Otherwise,
  this is stored in
  [`annotGeometries`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md).

- nucSeg:

  Similar to `cellSeg`, but for nuclei rather than whole cell.

## Usage

``` r
colGeometry(x, type = 1L, sample_id = 1L, withDimnames = TRUE)

colGeometry(
  x,
  type = 1L,
  sample_id = 1L,
  withDimnames = TRUE,
  translate = TRUE
) <- value

colGeometries(x, withDimnames = TRUE)

colGeometries(x, withDimnames = TRUE, translate = TRUE) <- value

colGeometryNames(x)

colGeometryNames(x) <- value

spotPoly(x, sample_id = 1L, withDimnames = TRUE)

spotPoly(x, sample_id = 1L, withDimnames = TRUE, translate = TRUE) <- value

centroids(x, sample_id = 1L, withDimnames = TRUE)

centroids(x, sample_id = 1L, withDimnames = TRUE, translate = TRUE) <- value

ROIPoly(x, sample_id = 1L, withDimnames = TRUE)

ROIPoly(x, sample_id = 1L, withDimnames = TRUE, translate = TRUE) <- value

cellSeg(x, sample_id = 1L, withDimnames = TRUE)

cellSeg(x, sample_id = 1L, withDimnames = TRUE, translate = TRUE) <- value

nucSeg(x, sample_id = 1L, withDimnames = TRUE)

nucSeg(x, sample_id = 1L, withDimnames = TRUE, translate = TRUE) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- type:

  An integer specifying the index or string specifying the name of the
  \*Geometry to query or replace. If missing, then the first item in the
  \*Geometries will be returned or replaced.

- sample_id:

  Sample ID to get or set geometries.

- withDimnames:

  Logical. If `TRUE`, then the dimnames (colnames or rownames) of the
  gene count matrix should correspond to row names of the `sf` data
  frames of interest.

- translate:

  Logical. Only used if
  [`removeEmptySpace`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/removeEmptySpace.md)
  has been run of the SFE object. If that's the case, this argument
  indicates whether the new value to be assigned to the geometry is in
  the coordinates prior to removal of empty space so it should be
  translated to match the new coordinates after removing empty space.
  Default to `TRUE`.

- value:

  Value to set. For `dimGeometry`, must be a `sf` data frame with the
  same number of rows as size in the dimension of interest, or an
  ordinary data frame that can be converted to such a `sf` data frame
  (see
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md)).
  For `dimGeometries`, must be a list of such `sf` or ordinary data
  frames.

## See also

\[dimGeometries()\], \[rowGeometries()\]

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> Error while performing HEAD request.
#>    Proceeding without cache information.
cgs <- colGeometries(sfe)
spots <- spotPoly(sfe)
```
